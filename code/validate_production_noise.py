"""Validate selected production spectral fields, preserving independent-run estimates.

Run from any directory. The protocol and seeds are fixed below; output goes to
data/analysis/revision_20261026/noise_validation.json. This does not run growth
simulations, regenerate archived fields, or fit an effective RG exponent.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time

os.environ["NUMBA_NUM_THREADS"] = "2"

import numpy as np

from bd_core import _amp_any, ffgn_weights, generate_spectral_eta
from check_generator_amplitude import ar_sum_covariance


ROOT = Path(__file__).resolve().parents[1]
LAGS = np.array([1, 2, 4, 8, 16, 32, 64, 128])
SPATIAL_LAGS = np.array([1, 2, 4, 8])
T = 65536
N_RUNS = 8
# Predeclared engineering tolerances, not simultaneous hypothesis-test p values.
SE_MULTIPLIER = 5.0
ABS_TOLERANCE = 0.003


def run_correlations(eta):
    """Circular products of population-centred signs, one estimate per run.

    The Gaussian source has population mean zero, so s=2*eta-1 has population
    mean zero and variance one. No estimated within-run mean is subtracted.
    Circular temporal products match the finite-period spectral target.
    """
    temporal, spatial = [], []
    for run in eta:
        sign = 2 * run.astype(np.int8) - 1
        temporal.append([float(np.mean(sign * np.roll(sign, -int(lag), axis=0)))
                         for lag in LAGS])
        spatial.append([float(np.mean(sign * np.roll(sign, -int(lag), axis=1)))
                        for lag in SPATIAL_LAGS])
    return np.array(temporal), np.array(spatial)


def summarise(temporal, spatial):
    n_runs = len(temporal)
    return {
        "C_time": temporal.mean(axis=0).tolist(),
        "C_time_err": (temporal.std(axis=0, ddof=1) / np.sqrt(n_runs)).tolist(),
        "C_space": spatial.mean(axis=0).tolist(),
        "C_space_err": (spatial.std(axis=0, ddof=1) / np.sqrt(n_runs)).tolist(),
        "per_run_C_time": temporal.tolist(),
        "per_run_C_space": spatial.tolist(),
    }


def finite_period_target(theta, spectrum):
    amplitude = _amp_any(theta, T, spectrum=spectrum)
    covariance = np.fft.irfft(amplitude ** 2, n=T)
    rho = covariance[LAGS] / covariance[0]
    return (2.0 / np.pi) * np.arcsin(np.clip(rho, -1.0, 1.0))


def checks(result, target):
    time_residual = np.abs(np.array(result["C_time"]) - target)
    time_bound = SE_MULTIPLIER * np.array(result["C_time_err"]) + ABS_TOLERANCE
    space_bound = SE_MULTIPLIER * np.array(result["C_space_err"]) + ABS_TOLERANCE
    return {
        "time_matches_finite_period_target": bool(np.all(time_residual <= time_bound)),
        "space_matches_zero": bool(np.all(np.abs(result["C_space"]) <= space_bound)),
        "max_abs_time_residual": float(time_residual.max()),
        "max_abs_space_correlation": float(np.abs(result["C_space"]).max()),
    }


def main():
    started = time.monotonic()
    core_path = ROOT / "code" / "bd_core.py"
    core_hash = hashlib.sha256(core_path.read_bytes()).hexdigest()
    series, legacy = [], None
    for spectrum in ("powerlaw", "fgn"):
        for theta in (0.25, 0.4):
            for L in (32, 64):
                seed = 2026102600 + len(series)
                temporal, spatial, old_temporal, old_spatial = [], [], [], []
                for run0, eta in generate_spectral_eta(
                        theta, L, T, N_RUNS, seed, spectrum=spectrum,
                        rows_per_chunk=2, max_bytes=2 * T * L):
                    assert run0 == len(temporal)
                    assert eta.shape[1:] == (T, L) and eta.dtype == np.uint8
                    time_values, space_values = run_correlations(eta)
                    temporal.extend(time_values)
                    spatial.extend(space_values)
                    if spectrum == "fgn" and theta == 0.25 and L == 32:
                        old = eta.transpose(0, 2, 1).reshape(eta.shape)
                        time_values, space_values = run_correlations(old)
                        old_temporal.extend(time_values)
                        old_spatial.extend(space_values)
                result = {"generator": spectrum, "theta": theta, "L": L, "T": T,
                          "n_runs": N_RUNS, "seed": seed,
                          **summarise(np.array(temporal), np.array(spatial))}
                target = finite_period_target(theta, spectrum)
                result["target_time"] = target.tolist()
                result["checks"] = checks(result, target)
                assert result["checks"]["time_matches_finite_period_target"], result
                assert result["checks"]["space_matches_zero"], result
                series.append(result)
                print(f"{spectrum} theta={theta} L={L}: "
                      f"C_time(1)={result['C_time'][0]:.6f}, "
                      f"target={target[0]:.6f}, "
                      f"max|C_space|={result['checks']['max_abs_space_correlation']:.6f}",
                      flush=True)
                if old_temporal:
                    legacy = {"theta": theta, "L": L, "T": T, "n_runs": N_RUNS,
                              "seed": seed, "generator": spectrum,
                              "operation": "correct_eta.transpose(0,2,1).reshape(correct_eta.shape)",
                              **summarise(np.array(old_temporal), np.array(old_spatial))}
                    legacy["checks"] = checks(legacy, target)
                    assert not legacy["checks"]["time_matches_finite_period_target"]
                    assert not legacy["checks"]["space_matches_zero"]

    weights, rates = ffgn_weights(0.0, 48)
    shared = ar_sum_covariance(weights, rates, np.r_[0, LAGS])
    independent = ar_sum_covariance(weights, rates, np.r_[0, LAGS], shared=False)
    infinite_fgn = {}
    for theta in (0.25, 0.4):
        H = theta + 0.5
        tau = LAGS.astype(float)
        rho = 0.5 * ((tau + 1) ** (2 * H) - 2 * tau ** (2 * H)
                     + (tau - 1) ** (2 * H))
        infinite_fgn[str(theta)] = ((2 / np.pi) * np.arcsin(rho)).tolist()
    final_core_hash = hashlib.sha256(core_path.read_bytes()).hexdigest()
    assert final_core_hash == core_hash, "bd_core.py changed during validation"
    output = {
        "lags": LAGS.tolist(), "spatial_lags": SPATIAL_LAGS.tolist(),
        "series": series, "legacy": legacy,
        "ar_zero": {"lags": LAGS.tolist(),
                    "C_shared": (shared[1:] / shared[0]).tolist(),
                    "C_independent": (independent[1:] / independent[0]).tolist()},
        "metadata": {
            "protocol": "production-noise-validation-v1",
            "noise_layout": "spectral-run-time-site-v2",
            "correlation_estimator": "per-run circular mean of (2*eta-1) products; population mean zero, variance one",
            "error_estimator": "sample standard deviation of 8 independent-run estimates, ddof=1, divided by sqrt(8)",
            "target_definition": "binary rho=(2/pi)*arcsin(irfft(_amp_any(theta,T,spectrum)**2)[lag]/irfft(_amp_any(theta,T,spectrum)**2)[0])",
            "target_boundary": "finite-period Gaussian spectrum, DC removed; fgn alias sum is truncated as in the production API; not the infinite-fGn covariance",
            "infinite_fgn_binary_reference_by_theta": infinite_fgn,
            "legacy_boundary": "this repository's old axis-layout fault only; same corrected fields as fgn theta=0.25 L=32 series",
            "ar_zero_boundary": "exact stationary raw-field normalised covariance of uniform-driven AR sums, theta=0 N=48 a=6 B=2; no Gaussian arcsine conversion; not the simulator's iid theta=0 shortcut",
            "checks": {"SE_multiplier": SE_MULTIPLIER, "absolute_tolerance": ABS_TOLERANCE,
                       "interpretation": "predeclared numerical smoke-check tolerances, not simultaneous confidence levels"},
            "rows_per_chunk": 2, "runs_per_chunk": 2, "NUMBA_NUM_THREADS": 2,
            "source_sha256": {"code/bd_core.py": core_hash,
                              "code/check_generator_amplitude.py": hashlib.sha256(
                                  (ROOT / "code" / "check_generator_amplitude.py").read_bytes()).hexdigest(),
                              "code/validate_production_noise.py": hashlib.sha256(
                                  Path(__file__).read_bytes()).hexdigest()},
            "elapsed_seconds": time.monotonic() - started,
        },
    }
    destination = ROOT / "data" / "analysis" / "revision_20261026" / "noise_validation.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {destination}", flush=True)


if __name__ == "__main__":
    main()
