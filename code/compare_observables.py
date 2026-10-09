"""
Compare candidate distributional (Hellinger) observables on the SAME data, to
establish which one is actually graded in theta and has an acceptable null
floor.

Candidates
  incr    : exact pmf of the nearest-neighbour increment, support {-1,0,+1}
            (binning-free)
  terr    : exact pmf of flat-terrace lengths, clipped at 48 (binning-free)
  diff    : histogram of standardised height differences at separation r
            (reported for comparison; suffers from support mismatch)
  height  : histogram of standardised heights (blind control)

Usage: python code/compare_observables.py --tag spec --L 128 256
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
from diagnostics import (hellinger_excess, hellinger_null_curve,
                         increment_pmf, terrace_length_pmf,
                         discrete_hellinger_excess)

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="spec")
ap.add_argument("--rule", default="bd")
ap.add_argument("--L", type=int, nargs="+", default=[128, 256])
ap.add_argument("--splits", type=int, default=80)
ap.add_argument("--r-frac", type=float, default=0.125)
args = ap.parse_args()

files = discover(args.tag, args.rule)
thetas_present = sorted({k[0] for k in files})

for L in args.L:
    if (0.0, L) not in files:
        print(f"L={L}: no reference")
        continue
    ref = load_result(files[(0.0, L)])
    ref_last = ref["profiles"][:, -1, :]
    r = max(2, int(args.r_frac * L))

    ref_inc = increment_pmf(ref_last)
    ref_terr = terrace_length_pmf(ref_last)

    # null floors
    hn_i = discrete_hellinger_excess(ref_inc, ref_inc, args.splits, 5)
    hn_t = discrete_hellinger_excess(ref_terr, ref_terr, args.splits, 5)
    hn_d, hsd_d = hellinger_null_curve(ref_last, 64, args.splits, 7, "diff", r)
    hn_h, hsd_h = hellinger_null_curve(ref_last, 64, args.splits, 7, "height", 1)
    print(f"=== L={L} (r={r})  null floors: incr={hn_i['Hnull_mean']:.5f}"
          f"+-{hn_i['Hnull_std']:.5f}  terr={hn_t['Hnull_mean']:.5f}"
          f"+-{hn_t['Hnull_std']:.5f}  diff={hn_d:.4f}  height={hn_h:.4f}")
    print(f"{'theta':>9} | {'incr dH':>9} {'z':>7} | {'terr dH':>9} {'z':>7} "
          f"| {'diff dH':>9} | {'height dH':>9} | {'p_flat':>7} {'<terr>':>7}")

    for th in [t for t in thetas_present if (t, L) in files]:
        cell = load_result(files[(th, L)])
        last = cell["profiles"][:, -1, :]
        ri = discrete_hellinger_excess(ref_inc, increment_pmf(last),
                                       args.splits, 5)
        rt = discrete_hellinger_excess(ref_terr, terrace_length_pmf(last),
                                       args.splits, 5)
        rd = hellinger_excess(ref_last, last, 64, args.splits, 7, None,
                              "diff", r)
        rh = hellinger_excess(ref_last, last, 64, args.splits, 7, None,
                              "height", 1)
        d = last.astype(np.int64)
        dd = d - np.roll(d, -1, axis=-1)
        pf = (dd == 0).mean()
        print(f"{th:9.4f} | {ri['dH']:9.5f} {ri['z']:7.2f} | {rt['dH']:9.5f} "
              f"{rt['z']:7.2f} | {rd['dH']:9.4f} | {rh['dH']:9.4f} | "
              f"{pf:7.4f} {np.mean([terrace_length_pmf(last)[i].argmax() for i in range(2)]):7.1f}")
    print()
