"""Fit binary, rather than latent, covariance and independently validate it.

Fit only the first 32 archived uniform-shared noise runs. The positive mixture
of circular Gaussian AR spectra defines an approximate binary match at sampled
lags, not an exact identity of complete processes. Independent uniform and
Gaussian validation streams use dense lags 0..64 plus the production long lags.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import fft, optimize
from numba import set_num_threads

from bd_core import ffgn_weights
from control_archives import atomic_savez
from innovation_controls import (ROOT, grow, noise_diagnostics, noise_field,
                                 gaussian_sign_covariance)


def spectral_basis(period):
    rates = np.r_[0., np.exp(-np.geomspace(6, 1/period, 60))]
    phase = np.exp(-2j*np.pi*np.fft.rfftfreq(period))
    spectra = 1/abs(1-rates[:, None]*phase)**2
    spectra[:, 0] = 0
    variances = fft.irfft(spectra, n=period, axis=-1)[:, 0]
    return spectra/variances[:, None]


def fit_binary_spectrum(target, lags, period):
    target, lags = np.asarray(target), np.asarray(lags)
    if target.shape != lags.shape or np.any(abs(target) > 1+1e-12) or lags[0] != 0:
        raise ValueError('need valid sign covariance at nonnegative lags including zero')
    basis = spectral_basis(period)
    cov = fft.irfft(basis, n=period, axis=-1)[:, lags].T
    weight = np.ones(len(lags)); weight[0] = 100.
    coefficients, _ = optimize.nnls(cov*weight[:, None],
                                     np.sin(np.pi/2*target)*weight, maxiter=10000)
    coefficients /= coefficients.sum()
    spectrum = coefficients @ basis
    predicted = gaussian_sign_covariance(spectrum, period)[lags]
    return spectrum, dict(coefficients=coefficients.tolist(),
                          max_training_residual=float(max(abs(predicted-target))),
                          predicted_training_covariance=predicted.tolist())


def sample_run(run, spec, amp, validation_lags):
    set_num_threads(1)
    seed = np.random.SeedSequence([20261101, int(round(spec['theta']*1e6)), spec['L'], run])
    gaussian_seed, uniform_seed = [int(s.generate_state(1)[0]) for s in seed.spawn(2)]
    rng = np.random.default_rng(gaussian_seed)
    period = spec['record_length']
    eta = np.empty((period, spec['L']), dtype=np.uint8)
    for start in range(0, spec['L'], 8):
        n = min(8, spec['L']-start)
        x = fft.irfft(fft.rfft(rng.standard_normal((n, period)), axis=-1)*amp,
                      n=period, axis=-1)
        eta[:, start:start+n] = (x > 0).T
    result = {**grow(eta, spec['T']), **noise_diagnostics(eta, validation_lags)}
    del eta
    w, r = ffgn_weights(spec['theta'], spec['N'])
    reference = noise_field(w, r, period, 16, uniform_seed, 'uniform_shared')
    diag = noise_diagnostics(reference, validation_lags)
    result['validation_reference_covariance'] = diag['sign_covariance']
    result['validation_reference_persistence'] = diag['persistence']
    return result


def run_binary_cell(source, outdir, workers=8):
    with np.load(source, allow_pickle=False) as z:
        spec = json.loads(str(z['spec']))
        training = z['sign_covariance'][:32, 3].mean(axis=0)
    period = spec['record_length']
    lags = np.unique(np.r_[np.arange(65), spec['lags']])
    spectrum, fit = fit_binary_spectrum(training, spec['lags'], period)
    metadata = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    driver_sha256=hashlib.sha256((ROOT/'code/innovation_controls.py').read_bytes()).hexdigest(),
                    source_file=source.relative_to(ROOT).as_posix(), training_runs=32,
                    validation='fresh independent uniform and Gaussian streams; 16 sites per run',
                    covariance_equivalence_margin=.005, fit=fit, spec=spec)
    path = outdir/source.name
    if path.exists():
        with np.load(path, allow_pickle=False) as z:
            if json.loads(str(z['metadata'])) != metadata:
                raise ValueError(f'incompatible cached binary cell: {path}')
        print(f'Validated cached binary {path.name}', flush=True)
        return
    t0 = time.time()
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, row in enumerate(pool.map(lambda run: sample_run(run, spec, np.sqrt(spectrum), lags),
                                          range(spec['runs']))):
            rows.append(row)
            if (i+1) % 16 == 0:
                print(f'Binary theta={spec["theta"]:.2f} L={spec["L"]} {i+1}/{spec["runs"]}, {time.time()-t0:.1f}s', flush=True)
    data = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    expected = gaussian_sign_covariance(spectrum, period)[lags]
    atomic_savez(path, metadata=json.dumps(metadata), validation_lags=lags,
                 expected_sign_covariance=expected, fitted_spectrum=spectrum,
                 **data, wall_seconds=time.time()-t0)
    print(f'Saved binary {path.name}', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, default=ROOT/'data/raw/innovation_20261031_v3')
    ap.add_argument('--out', type=Path, default=ROOT/'data/raw/innovation_binary_20261031_v3')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--wait', action='store_true', help='wait for each formal source cell')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for sites in [32, 64, 128, 256, 512]:
        for theta in [.1, .25, .4]:
            source = args.source/f'theta{theta:.2f}_L{sites}.npz'
            while args.wait and not source.exists():
                time.sleep(5)
            if source.exists():
                run_binary_cell(source, args.out, args.workers)


if __name__ == '__main__':
    main()
