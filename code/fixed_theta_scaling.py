"""
DECISIVE TEST: at fixed theta, does the deviation from the KPZ reference GROW
with system size?

This is the clean discriminator between the two scenarios, and it avoids the
statistical trap of the theta*(L) drift analysis.

If the model has a crossover length xi_x(theta) that grows as theta -> 0+, then
for a FIXED theta the surface is KPZ-like at scales below xi_x and long-range
above it.  Increasing L therefore resolves MORE of the long-range behaviour, so
the deviation from the theta = 0 reference must INCREASE with L:
        dH(theta, L) increasing in L   ->  genuine crossover (or threshold)
If instead the deviation is a purely local, one-point effect with no
scale-dependence, it is L-independent:
        dH(theta, L) constant in L     ->  no crossover: the observable merely
                                           tracks a smooth local parameter

Contrast with the theta*(L) drift: the split-half null floor shrinks with L
(more columns per run), so a fixed-significance criterion necessarily moves to
smaller theta as L grows EVEN IF the physics is scale-free.  The fixed-theta
test below is immune to that.

Usage: python code/fixed_theta_scaling.py
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="spec",
                help="analysis tag; its data/analysis/<tag>_bd/diagnostics.csv "
                     "must exist (run code/analyze.py first)")
ap.add_argument("--rule", default="bd")
args = ap.parse_args()


def load(tag, rule="bd"):
    p = os.path.join(ROOT, "data", "analysis", f"{tag}_{rule}", "diagnostics.csv")
    with open(p, newline="") as f:
        r = list(csv.DictReader(f))
    out = {}
    for k in r[0]:
        try:
            out[k] = np.array([float(x[k]) for x in r])
        except (ValueError, KeyError):
            pass
    return out


def main():
    d = load(args.tag, args.rule)
    Ls = sorted(set(d["L"].astype(int)))
    thetas = sorted(set(d["theta"]))
    print("=" * 104)
    print(f"dH(theta, L) at FIXED theta  [tag={args.tag}]")
    print("=" * 104)
    hdr = f"{'theta':>9} | " + " ".join(f"{'L='+str(L):>10}" for L in Ls) + \
          f" | {'ratio 512/32':>13}"
    print(hdr)
    print("-" * 104)
    for th in thetas:
        vals = []
        for L in Ls:
            m = (d["L"].astype(int) == L) & (np.abs(d["theta"] - th) < 1e-9)
            vals.append(d["dH"][m][0] if m.sum() else np.nan)
        v = np.array(vals)
        ok = np.isfinite(v)
        ratio = (v[ok][-1] / v[ok][0]) if ok.sum() >= 2 and abs(v[ok][0]) > 1e-9 else np.nan
        print(f"{th:9.4f} | " + " ".join(
            f"{x:10.5f}" if np.isfinite(x) else f"{'--':>10}" for x in v)
            + f" | {ratio:13.2f}" if np.isfinite(ratio) else
            f"{th:9.4f} | " + " ".join(
                f"{x:10.5f}" if np.isfinite(x) else f"{'--':>10}" for x in v) + f" | {'--':>13}")

    print("\n" + "=" * 104)
    print("Same, for the null floor (to show it shrinks with L, which is what")
    print("contaminates a fixed-significance onset drift)")
    print("=" * 104)
    print(f"{'L':>6} {'null floor':>12} {'null std':>10}")
    for L in Ls:
        m = d["L"].astype(int) == L
        print(f"{L:6d} {d['Hnull'][m][0]:12.5f} {d['Hnull_std'][m][0]:10.5f}")

    print("\n" + "=" * 104)
    print("RATIO TEST  (fixed theta): IQR-free summary of the L-dependence")
    print("=" * 104)
    print(f"{'theta':>9} {'dH(L=32)':>10} {'dH(L=512)':>10} {'trend':>12} {'significant?':>14}")
    for th in thetas:
        v, err, LL = [], [], []
        for L in Ls:
            m = (d["L"].astype(int) == L) & (np.abs(d["theta"] - th) < 1e-9)
            if m.sum():
                v.append(d["dH"][m][0]); err.append(d["dH_err"][m][0]); LL.append(L)
        v = np.array(v); err = np.array(err); LL = np.array(LL)
        if v.size >= 3:
            sl = np.polyfit(np.log(LL), v, 1)[0]     # d dH / d ln L
            errm = np.mean(err)
            verdict = "yes" if abs(sl) * np.log(LL[-1] / LL[0]) > 2 * errm else "no"
            trend = "grows" if sl > 0 else "shrinks"
            print(f"{th:9.4f} {v[0]:10.5f} {v[-1]:10.5f} {trend:>12} {verdict:>14}"
                  f"   (d dH/d lnL = {sl:+.5f})")


if __name__ == "__main__":
    main()
