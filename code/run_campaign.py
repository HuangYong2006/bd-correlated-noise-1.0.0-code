"""
Campaign runner (literature-standard correlated BD).

Parallelism is provided by numba ``prange`` inside ``bd_core.simulate``
(threads), NOT by multiprocessing: the DSH file sandbox blocks the named pipes
that joblib/loky require, so process-based parallelism is unavailable here.

The number of sweeps is scaled with the system size as  T = t_mult * L^1.5 ,
which is proportional to the measured KPZ saturation time t_x ~ 0.6 L^1.5.

Usage
-----
python code/run_campaign.py --tag main --L 128 256 512 1024 2048 --t-mult 20 \
    --runs 160 --thetas 0 0.05 0.10 0.15 0.1667 0.20 0.25 0.30 0.40
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from bd_core import (SimSpec, simulate, simulate_spectral, save_result,
                     RULE_BD, RULE_RSOS, smooth_ceil)

DATA = os.path.join(ROOT, "data", "raw")


def cell_file(tag, rule_name, theta, L):
    th = f"{theta:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
    return os.path.join(DATA, tag, f"{rule_name}_{th}_L{L}.npz")


def theta_str(theta):
    return f"{theta:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--thetas", type=float, nargs="+", required=True)
    ap.add_argument("--L", type=int, nargs="+", required=True)
    ap.add_argument("--t-mult", type=float, default=20.0,
                    help="T = t_mult * L^1.5 (sweeps)")
    ap.add_argument("--t-cap", type=int, default=8_000_000)
    ap.add_argument("--runs", type=int, default=160)
    ap.add_argument("--rule", choices=["bd", "rsos"], default="bd")
    ap.add_argument("--generator", choices=["ffgn", "spectral"], default="ffgn",
                    help="ffgn: Mandelbrot FFGN, the published-literature "
                         "generator (discontinuous at theta=0); "
                         "spectral: exact power-law spectral filtering, "
                         "continuous at theta=0 (correct for a theta-scan)")
    ap.add_argument("--spectrum", choices=["powerlaw", "fgn"], default="powerlaw",
                    help="filter amplitude used when --generator spectral.  "
                         "'powerlaw' is |f|^(-theta), the family used for the "
                         "'spec' campaign; its binarised theta axis is accurate "
                         "to ~10%% and it under-correlates at large theta.  "
                         "'fgn' is the exact fractional-Gaussian-noise spectral "
                         "density, whose binarised theta axis is accurate to "
                         "~1%% and which is exactly white at theta = 0.")
    ap.add_argument("--n-check", type=int, default=32)
    ap.add_argument("--n-w", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20240)
    ap.add_argument("--skip-existing", action="store_true")
    args = ap.parse_args()

    rule = RULE_BD if args.rule == "bd" else RULE_RSOS
    outdir = os.path.join(DATA, args.tag)
    os.makedirs(outdir, exist_ok=True)

    tasks = [(th, L) for L in args.L for th in args.thetas]
    print(f"{len(tasks)} cells: {len(args.thetas)} thetas x {len(args.L)} sizes x "
          f"{args.runs} runs   (T = {args.t_mult} * L^1.5, cap {args.t_cap})", flush=True)

    t_all = time.time()
    manifest = {"tag": args.tag, "rule": args.rule, "runs": args.runs,
                "t_mult": args.t_mult, "cells": []}

    for (th, L) in tasks:
        T = int(min(args.t_cap, round(args.t_mult * L ** 1.5)))
        if args.generator == "spectral":
            T = smooth_ceil(T)          # keep the FFT length 5-smooth
        path = cell_file(args.tag, args.rule, th, L)
        if args.skip_existing and os.path.exists(path):
            print(f"  skip {os.path.basename(path)}", flush=True)
            continue
        seed = int(np.random.SeedSequence(
            [args.seed, int(round(th * 1e6)) + 10 ** 7, L]).generate_state(1)[0] % (2 ** 31))
        spec = SimSpec(theta=th, L=L, T=T, n_runs=args.runs, rule=rule,
                       n_check=args.n_check, n_w=args.n_w, seed=seed,
                       spectrum=args.spectrum)
        t0 = time.time()
        res = (simulate_spectral(spec) if args.generator == "spectral"
               else simulate(spec))
        save_result(path, res)
        manifest["cells"].append({"theta": th, "L": L, "T": T, "n_runs": args.runs,
                                  "file": os.path.basename(path),
                                  "wall": round(time.time() - t0, 2)})
        print(f"  [{len(manifest['cells']):3d}/{len(tasks)}] theta={th:+.4f} L={L:5d} "
              f"T={T:9d}  {time.time()-t0:7.1f}s  -> {os.path.basename(path)}",
              flush=True)

    with open(os.path.join(outdir, f"manifest_{args.rule}.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nTOTAL: {time.time()-t_all:.1f}s", flush=True)


if __name__ == "__main__":
    main()
