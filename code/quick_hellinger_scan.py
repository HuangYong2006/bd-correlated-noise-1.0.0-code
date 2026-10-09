"""
Quick targeted check: is the distributional (Hellinger) diagnostic GRADED in
theta once the generator is theta-continuous?

With the FFGN generator the observable jumps to its near-maximum at the first
nonzero theta on the grid and is then flat, because FFGN(theta->0) is already
strongly correlated.  A theta-continuous generator should instead produce a
smooth rise, which is what a finite-size-scaling drift analysis needs.

This script computes only the Hellinger diagnostic (no Isolation Forest), so it
is cheap and can be run while a campaign is still filling in cells.

Usage: python code/quick_hellinger_scan.py --tag spec --L 64 128 256
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analyze import discover
from bd_core import load_result
from diagnostics import hellinger_excess, hellinger_null_curve

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="spec")
ap.add_argument("--rule", default="bd")
ap.add_argument("--L", type=int, nargs="+", default=[64, 128, 256])
ap.add_argument("--bins", type=int, default=64)
ap.add_argument("--splits", type=int, default=80)
ap.add_argument("--r-frac", type=float, default=0.125)
args = ap.parse_args()

files = discover(args.tag, args.rule)
thetas_present = sorted({k[0] for k in files})
print(f"tag={args.tag}  {len(files)} cells  thetas={thetas_present}\n")

for L in args.L:
    if (0.0, L) not in files:
        print(f"L={L}: no theta=0 reference yet, skipping")
        continue
    ref = load_result(files[(0.0, L)])
    ref_last = ref["profiles"][:, -1, :]
    r = max(2, int(args.r_frac * L))
    hn, hsd = hellinger_null_curve(ref_last, args.bins, args.splits, 7, "diff", r)
    print(f"=== L={L}  probe r={r}  null floor = {hn:.4f} +- {hsd:.4f} ===")
    print(f"{'theta':>9} {'dH':>9} {'z':>8} {'p_flat':>8}")
    ths = [t for t in thetas_present if (t, L) in files]
    for th in ths:
        cell = load_result(files[(th, L)])
        last = cell["profiles"][:, -1, :]
        res = hellinger_excess(ref_last, last, args.bins, args.splits, 7, None,
                               "diff", r)
        d = last.astype(np.int64)
        dd = d - np.roll(d, -1, axis=-1)
        print(f"{th:9.4f} {res['dH']:9.4f} {res['z']:8.2f} {(dd == 0).mean():8.4f}")
    print()
