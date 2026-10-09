"""Determine how many realisations the diagnostics actually used per cell.

The review found an apparent contradiction: the text says run counts were
"matched at 96", while the stored cells hold 192 (L<=256) and 96 (L=512) for the
spectral campaigns and 256 for the FFGN campaign.  There are two possibilities:
  (a) 96 is the real number of stored realisations for every corrected cell; or
  (b) 96 is a common subsample applied by analyze.py (--if-runs 96 --hell-runs 96).

The split-half null floor depends on n_half = min(n_ref, n_th)//2, so the value
recorded in diagnostics.csv identifies the sample size that was used.  This
script recomputes the null floor for several candidate sample sizes and compares.
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from diagnostics import discrete_hellinger_excess, increment_pmf  # noqa: E402

L = 256


def csv_row(tag, theta, L):
    p = os.path.join(ROOT, "data", "analysis", f"{tag}_bd", "diagnostics.csv")
    with open(p, newline="") as f:
        for r in csv.DictReader(f):
            if abs(float(r["theta"]) - theta) < 1e-9 and int(float(r["L"])) == L:
                return r
    return None


def profiles(tag, theta):
    tok = "p0d%04d" % round(theta * 1e4)
    p = os.path.join(ROOT, "data", "raw", tag, f"bd_{tok}_L{L}.npz")
    return np.load(p, allow_pickle=False)["profiles"]


for tag, stored in (("fgnspec", 192), ("spec", 192), ("main", 256)):
    row = csv_row(tag, 0.0, L)
    if row is None:
        print(f"{tag}: no theta=0 row at L={L}")
        continue
    recorded = float(row["Hnull"])
    print(f"\n=== {tag}  (theta=0, L={L})  stored realisations = {stored}")
    print(f"    recorded null floor in diagnostics.csv = {recorded:.6f}")
    prof = profiles(tag, 0.0)
    last = prof[:, -1, :]
    full = last.shape[0]
    for k in (96, 128, 192, 256):
        if k > full:
            continue
        pmf = increment_pmf(last[:k])
        res = discrete_hellinger_excess(pmf, pmf, 120, 0)
        flag = "  <-- matches the recorded value" if abs(
            res["Hnull_mean"] - recorded) < 2e-4 else ""
        print(f"      k={k:4d} runs (n_half={k//2:3d}) -> null floor "
              f"{res['Hnull_mean']:.6f}{flag}")
