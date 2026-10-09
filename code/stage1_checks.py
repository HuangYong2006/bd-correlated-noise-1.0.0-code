"""
Stage-1 checks.

A. Generator validation: does the quantile-mapped spectral-filtering sequence
   really have C(tau) ~ tau^(2*theta-1), and is its marginal exactly uniform?
B. Speed benchmark: how long does one BD ensemble take as a function of L, T?
C. KPZ sanity check at theta = 0: do we recover beta = 1/3, alpha = 1/2, z = 3/2?

Run:  python code/stage1_checks.py
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bd_core import (SimSpec, generate_columns, simulate, empirical_correlation,
                     RULE_BD)

rng = np.random.default_rng(0)


def logbin(x, y, nbins=24):
    """Log-spaced binning average of y over x>0."""
    m = x > 0
    x, y = x[m], y[m]
    edges = np.geomspace(x.min(), x.max(), nbins + 1)
    idx = np.digitize(x, edges) - 1
    xs, ys, ns = [], [], []
    for k in range(nbins):
        s = idx == k
        if s.sum() > 0:
            xs.append(x[s].mean())
            ys.append(y[s].mean())
            ns.append(int(s.sum()))
    return np.array(xs), np.array(ys), np.array(ns)


def fit_slope(x, y, lo, hi):
    m = (x >= lo) & (x <= hi) & np.isfinite(y) & (y > 0)
    if m.sum() < 5:
        return np.nan, np.nan
    p = np.polyfit(np.log(x[m]), np.log(y[m]), 1)
    return p[0], m.sum()


# --------------------------------------------------------------------------
print("=" * 78)
print("A. GENERATOR VALIDATION   (expected slope of C(tau) = 2*theta - 1)")
print("=" * 78)

T = 1 << 19
NSEQ = 24
max_lag = 20000
thetas = [0.0, 0.05, 0.10, 0.1666667, 0.20, 0.30, 0.40]

t0 = time.time()
# theta = 0 reference: i.i.d. uniform
u0 = rng.random((NSEQ, T))
lab = []
for th in thetas:
    if th == 0.0:
        ac = np.mean([empirical_correlation(u0[i], max_lag) for i in range(NSEQ)], axis=0)
    else:
        # For the correlation check we use the pre-quantisation Gaussian sequence
        # (this isolates the spectral recipe from the lattice quantisation).
        acl = []
        from bd_core import _spectral_amplitude
        from scipy.special import ndtr
        amp = _spectral_amplitude(th, T)
        ss = np.random.SeedSequence(2000 + int(th * 1e5))
        r = np.random.default_rng(ss)
        u = r.standard_normal((NSEQ, T))
        U = np.fft.rfft(u, axis=1)
        U *= amp
        s = np.fft.irfft(U, n=T, axis=1)
        s -= s.mean(axis=1, keepdims=True)
        s /= s.std(axis=1, keepdims=True)
        for i in range(NSEQ):
            acl.append(empirical_correlation(s[i], max_lag))
        ac = np.mean(acl, axis=0)
    xs, ys, ns = logbin(np.arange(1, max_lag), ac[1:])
    sl, npts = fit_slope(xs, np.abs(ys), 100, 8000)
    lab.append((th, 2 * th - 1, sl, ac[1], ac[10], ac[100]))
    print(f"  theta={th:+.4f}  expected={2*th-1:+.3f}  measured C(tau) slope={sl:+.3f}"
          f"   C(1)={ac[1]:+.4f} C(10)={ac[10]:+.4f} C(100)={ac[100]:+.5f}")

print(f"  [generator validation took {time.time()-t0:.1f}s]")

# marginal uniformity check
L_check = 64
cols_all = []
for i0, cols in generate_columns(0.2, L_check, 1 << 18, 4, seed=7):
    cols_all.append(cols.ravel())
c = np.concatenate(cols_all)
cnt = np.bincount(c, minlength=L_check)
chi2 = ((cnt - cnt.mean()) ** 2 / cnt.mean()).sum()
print(f"  marginal uniformity at theta=0.2, L={L_check}: "
      f"chi2/dof = {chi2/(L_check-1):.4f}  (expected ~1)")

# ---- correlation of the ACTUAL quantised column-index sequence used by BD ----
print()
print("  Quantised column-index sequence (the real physical noise), L=4096:")
Lq = 4096
for th in [0.0, 0.10, 0.1666667, 0.30]:
    acl = []
    n_seq = 8
    for i0, cols in generate_columns(th, Lq, 1 << 18, n_seq, seed=31):
        for i in range(cols.shape[0]):
            acl.append(empirical_correlation(cols[i].astype(np.float64), max_lag))
    ac = np.mean(acl, axis=0)
    xs, ys, ns = logbin(np.arange(1, max_lag), ac[1:])
    sl, npts = fit_slope(xs, np.abs(ys), 100, 8000)
    print(f"    theta={th:+.4f}  expected slope={2*th-1:+.3f}  measured={sl:+.3f}"
          f"   C(1)/C(10)/C(100) = {ac[1]:+.5f} / {ac[10]:+.5f} / {ac[100]:+.6f}")

# --------------------------------------------------------------------------
print()
print("=" * 78)
print("B. SPEED BENCHMARK")
print("=" * 78)

for (L, Tt, nr) in [(1024, 1_000_000, 16), (4096, 1_000_000, 8), (1024, 10_000_000, 4)]:
    for th in [0.0, 0.2]:
        spec = SimSpec(theta=th, L=L, T=Tt, n_runs=nr, n_w=4000)
        res = simulate(spec, verbose=False)
        per_run = res["wall_time"] / nr
        print(f"  L={L:6d} T={Tt:9d} theta={th:.1f} n={nr:3d}: "
              f"{res['wall_time']:7.2f}s total, {per_run*1e6/Tt:6.2f} ns/step/run")

# --------------------------------------------------------------------------
print()
print("=" * 78)
print("C. KPZ SANITY CHECK AT theta = 0   (expect beta=1/3, alpha=1/2)")
print("=" * 78)

for L in [512, 1024, 2048]:
    Tt = 1_000_000
    nr = 64
    spec = SimSpec(theta=0.0, L=L, T=Tt, n_runs=nr, n_check=24, n_w=200000)
    res = simulate(spec, verbose=False)
    wt, w2 = res["w_times"], res["w_series"]
    W = np.sqrt(np.mean(w2.astype(np.float64), axis=0))
    Wsat = W[-1]
    tsat = L ** 1.5
    slopes = []
    for lo, hi in [(30, 300), (100, 1000), (300, 3000), (1000, 10000), (3000, 30000)]:
        s, n = fit_slope(wt, W, lo, hi)
        slopes.append(f"[{lo},{hi}]:{s:.3f}")
    print(f"  L={L:5d}  W_sat={Wsat:8.3f}  t_sat(est)={tsat:9.0f}   beta " + "  ".join(slopes))

print()
print("  alpha check from W_sat vs L is done properly in the campaign analysis.")
print("DONE")
