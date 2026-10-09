"""
Direct diagnostic of the theta = 0 dynamics: W(t) over 6 decades for both the
true ballistic rule and the RSOS rule, for several L.  Purpose: verify that the
simulator reproduces the known KPZ exponents beta = 1/3 (growth) and
alpha = 1/2 (saturation, W_sat ~ L^alpha) before any production campaign.

Run:  python code/debug_bd.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bd_core import _bd_ensemble, RULE_BD, RULE_RSOS

T = 10_000_000
N_RUNS = 48
L_LIST = [128, 256, 512, 1024, 2048, 4096]
CHECK = np.unique(np.round(np.geomspace(10, T, 40)).astype(np.int64))


def run_zero_theta(L, rule, n_runs=N_RUNS, T=T, check=CHECK, seed=99):
    rng = np.random.default_rng(seed)
    cols = rng.integers(0, L, size=(n_runs, T)).astype(np.uint16)
    profiles = np.zeros((n_runs, len(check), L), dtype=np.int32)
    wseries = np.zeros((n_runs, 1), dtype=np.float32)
    _bd_ensemble(cols, L, rule, T + 1, check, profiles, wseries)
    p = profiles.astype(np.float64)
    W = p.std(axis=-1)                      # (n_runs, n_check)
    return check, W


def slope(t, w, lo, hi):
    m = (t >= lo) & (t <= hi)
    if m.sum() < 4:
        return np.nan
    return np.polyfit(np.log(t[m]), np.log(w[m]), 1)[0]


for rule, name in [(RULE_BD, "BD (max(h_im,h_i+1,h_ip))"),
                   (RULE_RSOS, "RSOS (max(h_im,h_ip)+1)")]:
    print("=" * 110)
    print(f"RULE: {name}")
    print("=" * 110)
    wsat = {}
    for L in L_LIST:
        t, W = run_zero_theta(L, rule)
        wm = W.mean(axis=0)
        wsat[L] = wm[-1]
        b1 = slope(t, wm, 100, 1000)
        b2 = slope(t, wm, 1000, 10000)
        b3 = slope(t, wm, 10 ** 4, 10 ** 5)
        b4 = slope(t, wm, 10 ** 5, 10 ** 6)
        print(f"  L={L:5d}  W(1e2)={wm[np.argmin(abs(t-100))]:8.3f}"
              f"  W(1e4)={wm[np.argmin(abs(t-10**4))]:8.3f}"
              f"  W(1e6)={wm[np.argmin(abs(t-10**6))]:8.3f}"
              f"  W(1e7)={wm[-1]:8.3f}   beta: {b1:.3f} {b2:.3f} {b3:.3f} {b4:.3f}")
    print()
    Ls = np.array(sorted(wsat))
    Ws = np.array([wsat[L] for L in Ls])
    print("  saturated W_sat vs L (fit over the largest 4 sizes):")
    p = np.polyfit(np.log(Ls[-4:]), np.log(Ws[-4:]), 1)
    print(f"    alpha_fit = {p[0]:.4f}   (expected 0.5)")
    for L, w in zip(Ls, Ws):
        print(f"      L={L:5d}  W_sat={w:8.3f}   W_sat/sqrt(L)={w/np.sqrt(L):.4f}")
    print()

print("DONE")
