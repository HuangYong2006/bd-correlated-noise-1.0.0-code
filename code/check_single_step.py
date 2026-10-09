"""
Structural check of the 0/1-noise BD rule.

Claim to verify: with a binarised (0/1) noise and a parallel update,
    h(x,t+1) = max[ h(x,t) + eta'(x,t), h(x-1,t), h(x+1,t) ]
starting from a flat initial condition, the single-step constraint
|h(x,t) - h(x+1,t)| <= 1 is preserved for all t, so the model has NO overhang
degrees of freedom: the height-increment distribution has support on {-1,0,+1}
only and is completely described by two numbers.

Argument: given |Delta h| <= 1 everywhere at time t,
  * h(x-1,t), h(x+1,t) <= h(x,t) + 1  =>  h(x,t+1) <= h(x,t) + 1
  * h(x,t+1) >= h(x,t) + eta' >= h(x,t)
so h(x,t+1) - h(x,t) in {0,1}. If a neighbour is higher, that neighbour's value
is taken, so a site can never fall behind by more than one either; both effects
together keep |Delta h| <= 1.

We verify this numerically on the production cells and also measure the two
numbers that fully characterise the morphology in this model:
  p_flat = frac{Delta h = 0},  p_up = frac{Delta h = +1},  p_down = frac{Delta h = -1}.

Usage: python code/check_single_step.py --tag main --rule bd
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from bd_core import load_result


def discover(tag, rule):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = f
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main")
    ap.add_argument("--rule", default="bd")
    ap.add_argument("--max-cells", type=int, default=12)
    args = ap.parse_args()

    files = discover(args.tag, args.rule)
    keys = sorted(files)[:args.max_cells]
    print(f"{'theta':>8} {'L':>6} {'max|dh|':>8} {'p_flat':>8} {'p_up':>8} "
          f"{'p_down':>8} {'skew_d':>8} {'note'}")
    print("-" * 78)
    worst = 0
    for k in keys:
        theta, L = k
        r = load_result(files[k])
        h = r["profiles"][:, -1, :].astype(np.int64)
        d = h - np.roll(h, -1, axis=-1)
        m = int(np.abs(d).max())
        worst = max(worst, m)
        flat = float((d == 0).mean())
        up = float((d == +1).mean())
        dn = float((d == -1).mean())
        sk = float(((d - d.mean()) ** 3).mean() / max(d.std() ** 3, 1e-30))
        note = "single-step OK" if m <= 1 else "!! VIOLATION"
        print(f"{theta:8.4f} {L:6d} {m:8d} {flat:8.4f} {up:8.4f} {dn:8.4f} {sk:8.3f} {note}")

    print("-" * 78)
    print(f"worst |Delta h| over inspected cells: {worst}")
    if worst <= 1:
        print("CONFIRMED: the 0/1-noise parallel-update BD rule is single-step;")
        print("the increment density has support only on {-1,0,+1} and is fully")
        print("specified by p_flat and the up/down asymmetry.  'Faceting' in this")
        print("model therefore cannot involve overhangs: it must appear as large")
        print("flat terraces separated by unit steps.")
    else:
        print("NOT CONFIRMED: overhangs with |Delta h| > 1 do occur.")


if __name__ == "__main__":
    main()
