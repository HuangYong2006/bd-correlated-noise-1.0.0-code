"""Run-level analysis of the regenerated, versioned manuscript campaigns.

The reported alpha is a finite-size slope of late-time RMS width. It is not an
asymptotic exponent or a critical-threshold estimate. Old spectral archives are
rejected before analysis. Random partitions calibrate Hellinger's finite-sample
floor; an outer independent-run bootstrap supplies its estimator uncertainty.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from bd_core import load_result
from diagnostics import discrete_hellinger_excess, increment_pmf

ROOT = Path(__file__).resolve().parents[1]
TAGS = {"shared": "rev20261026_shared", "powerlaw": "rev20261026_powerlaw",
        "fgn": "rev20261026_fgn"}
SIZES = [32, 64, 128, 256, 512]
THETAS = [0.0, 0.10, 1 / 6, 0.20, 0.25, 0.30, 0.40]


def late_width_samples(w_series, tail_fraction=0.2):
    """Return one mean squared width per run in adjacent final-time windows."""
    count = max(1, int(w_series.shape[1] * tail_fraction))
    return (w_series[:, -count:].mean(axis=1),
            w_series[:, -2 * count:-count].mean(axis=1))


def require_corrected_layout(result, generator):
    if generator != "shared" and result.get("noise_layout") != "spectral-run-time-site-v2":
        raise ValueError("spectral noise_layout must be spectral-run-time-site-v2")


def fit_size_exponent(squared_by_L, n_boot=500, seed=7):
    """Fit log(sqrt(mean(W^2))) vs log(L), resampling independent runs."""
    sizes = sorted(squared_by_L)
    x = np.log(sizes)
    centered = x - x.mean()
    rng = np.random.default_rng(seed)
    central = np.array([np.log(np.sqrt(np.mean(squared_by_L[L]))) for L in sizes])
    draws = np.empty((n_boot, len(sizes)))
    for k, L in enumerate(sizes):
        values = np.asarray(squared_by_L[L])
        indices = rng.integers(len(values), size=(n_boot, len(values)))
        draws[:, k] = np.log(np.sqrt(values[indices].mean(axis=1)))
    slopes = draws @ centered / (centered @ centered)
    pairs = []
    for k, (a, b) in enumerate(zip(sizes[:-1], sizes[1:])):
        scale = np.log(b / a)
        pairs.append({"L_low": a, "L_high": b, "midpoint": float(np.sqrt(a * b)),
                      "alpha": float((central[k + 1] - central[k]) / scale),
                      "alpha_err": float(((draws[:, k + 1] - draws[:, k]) / scale).std(ddof=1))})
    return {"alpha": float(central @ centered / (centered @ centered)),
            "alpha_err": float(slopes.std(ddof=1)), "alpha_pairs": pairs}


def read_campaign(tag, generator):
    cells = {}
    for path in sorted((ROOT / "data/raw" / tag).glob("bd_*.npz")):
        result = load_result(str(path))
        require_corrected_layout(result, generator)
        spec = result["spec"]
        if spec["n_runs"] != 64:
            raise ValueError(f"unmatched independent-run count in {path}")
        tail, previous = late_width_samples(result["w_series"].astype(np.float64))
        profile = result["profiles"][:, -1, :].astype(np.int64)
        increments = profile - np.roll(profile, -1, axis=1)
        if np.max(np.abs(increments)) > 1:
            raise ValueError(f"single-step invariant failed in {path}")
        key = (round(float(spec["theta"]), 4), int(spec["L"]))
        cells[key] = {"tail": tail, "previous": previous,
                      "pflat": (increments == 0).mean(axis=1),
                      "pmf": increment_pmf(profile), "spec": spec,
                      "raw_file": path.relative_to(ROOT).as_posix(),
                      "raw_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                      "noise_layout": result["noise_layout"]}
    expected = {(round(t, 4), L) for t in THETAS for L in SIZES}
    if set(cells) != expected:
        raise ValueError(f"{tag}: campaign incomplete or unexpected cells; missing {sorted(expected - set(cells))}")
    return cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--hell-bootstrap", type=int, default=200)
    parser.add_argument("--splits", type=int, default=120)
    parser.add_argument("--outdir", type=Path, default=ROOT / "data/analysis/revision_20261026")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows, exponents, sources = [], [], []
    for generator, tag in TAGS.items():
        data = read_campaign(tag, generator)
        for i, ((theta, L), cell) in enumerate(sorted(data.items())):
            ref = data[(0.0, L)]["pmf"]
            diag = discrete_hellinger_excess(ref, cell["pmf"], args.splits,
                                            seed=20261026 + i,
                                            n_boot=args.hell_bootstrap if theta != 0 else 0)
            rng = np.random.default_rng(20261026 + i)
            idx = rng.integers(64, size=(args.bootstrap, 64))
            width_samples = np.sqrt(cell["tail"][idx].mean(axis=1))
            p_samples = cell["pflat"][idx].mean(axis=1)
            ratio_samples = np.sqrt(cell["tail"][idx].mean(axis=1) /
                                    cell["previous"][idx].mean(axis=1))
            row = {"generator": generator, "theta": theta, "L": L,
                   "W": float(np.sqrt(cell["tail"].mean())),
                   "W_err": float(width_samples.std(ddof=1)),
                   "p_flat": float(cell["pflat"].mean()),
                   "p_flat_err": float(p_samples.std(ddof=1)),
                   "dH": float(diag["dH"]) if theta != 0 else 0.0,
                   "dH_err": float(diag["dH_err"]) if theta != 0 else 0.0,
                   "Hnull_mean": float(diag["Hnull_mean"]),
                   "Hnull_std": float(diag["Hnull_std"]),
                   "late_ratio": float(np.sqrt(cell["tail"].mean() / cell["previous"].mean())),
                   "late_ratio_err": float(ratio_samples.std(ddof=1)),
                   "T": cell["spec"]["T"], "n_runs": 64,
                   "raw_file": cell["raw_file"], "noise_layout": cell["noise_layout"]}
            rows.append(row)
            sources.append({"raw_file": cell["raw_file"], "sha256": cell["raw_sha256"],
                            "spec": cell["spec"], "noise_layout": cell["noise_layout"]})
            print(f"{generator:8} theta={theta:.4f} L={L:3} W={row['W']:.4f} pflat={row['p_flat']:.4f}", flush=True)
        for i, theta in enumerate(THETAS):
            samples = {L: data[(round(theta, 4), L)]["tail"] for L in SIZES}
            fit = fit_size_exponent(samples, args.bootstrap, seed=901 + i)
            exponents.append({"generator": generator, "theta": round(theta, 4), **fit})
    summary = {"protocol": {"tags": TAGS, "n_runs": 64, "Ls": SIZES, "thetas": THETAS,
                            "width_window": "last 20 percent of stored width times",
                            "previous_window": "preceding 20 percent",
                            "bootstrap": args.bootstrap, "hell_bootstrap": args.hell_bootstrap,
                            "splits": args.splits, "zero_theta_dH": "zero by reference definition",
                            "interpretation": "finite-size late-time slopes; no theta_c inference",
                            "forest": "not used for the manuscript's physical conclusions"},
               "exponents": exponents, "cells": rows, "sources": sources}
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    with (args.outdir / "observables.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    with (args.outdir / "exponents.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["generator", "theta", "alpha", "alpha_err"])
        writer.writeheader(); writer.writerows({k: row[k] for k in writer.fieldnames} for row in exponents)
    print(f"Wrote {len(rows)} verified cells and {len(exponents)} finite-size fits", flush=True)


if __name__ == "__main__":
    main()
