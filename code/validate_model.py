"""
Validation of the literature-standard model (parallel update, binarised
long-range temporally correlated noise).

A. Noise validation: does the FFGN recipe give <eta'(t)eta'(t+tau)> ~ tau^(2th-1)?
   How does it compare with exact power-law spectral filtering at the same
   nominal theta?  (This is the generator-dependence question that could
   explain the spread of published theta_c values.)
B. theta = 0 dynamics: KPZ exponents beta = 1/3, alpha = 1/2, and a measurement
   of the saturation time t_x(L) needed to plan the production campaign.
C. Timing.

Run:  python code/validate_model.py
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from bd_core import (SimSpec, simulate, ffgn_weights, ffgn_n_components,
                     fft_noise_field, _u01, RULE_BD, RULE_RSOS)

# --------------------------------------------------------------------------
# A. noise validation
# --------------------------------------------------------------------------

@njit(cache=True)
def ffgn_binary(n_sites, T, Wn, rn, seed):
    """Binarised FFGN noise, (T, n_sites) uint8, independent sites."""
    N = Wn.shape[0]
    out = np.empty((T, n_sites), dtype=np.uint8)
    for x in range(n_sites):
        X = np.zeros(N)
        base = np.uint64(seed + 1) * np.uint64(2654435761) + np.uint64(x) * np.uint64(0x9E3779B97F4A7C15)
        for t in range(T):
            z = _u01(base + np.uint64(t)) - 0.5
            s = 0.0
            for n in range(N):
                xn = rn[n] * X[n] + z
                X[n] = xn
                s += Wn[n] * xn
            out[t, x] = 1 if s > 0.0 else 0
    return out


def fft_autocorr(x):
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conjugate(F), n=nfft)[:n]
    return ac / ac[0]


def bin_corr(x, edges):
    """Correlation of a binary sequence via lag binning (fast, vectorised)."""
    return fft_autocorr(x.astype(np.float64))


def slope_of(tau, c, lo, hi):
    m = (tau >= lo) & (tau <= hi)
    return np.polyfit(np.log(tau[m]), np.log(np.abs(c[m])), 1)[0]


print("=" * 104)
print("A. NOISE GENERATOR VALIDATION")
print("=" * 104)

T_noise = 1 << 20
max_lag = 1 << 14
thetas = [0.05, 0.10, 0.15, 1 / 6, 0.20, 0.25, 0.30, 0.40]
N_SITES = 24

print(f"T = 2^{int(np.log2(T_noise))} sweeps, {N_SITES} sites averaged, "
      f"slope window tau in [64, 4096]")
print(f"\n{'theta':>7} {'2th-1':>7} | {'FFGN slope':>10} {'th_eff':>7} {'P(1)':>6} "
      f"| {'FFT slope':>10} {'th_eff':>7} | {'d(th_eff)':>9}")
print("-" * 104)

noise_rows = []
for th in thetas:
    N = ffgn_n_components(T_noise)
    Wn, rn = ffgn_weights(th, N)
    eta = ffgn_binary(N_SITES, T_noise, Wn, rn, seed=101 + int(th * 1000))
    ac_f = np.mean([bin_corr(eta[:, x], None) for x in range(N_SITES)], axis=0)
    tau = np.arange(1, max_lag)
    s_f = slope_of(tau, ac_f[1:max_lag], 64, 4096)
    th_f = (1.0 + s_f) / 2.0
    p1 = eta.mean()

    # exact spectral filtering at reduced T (memory/cost)
    T_fft = 1 << 17
    eta2 = fft_noise_field(th, 8, T_fft, seed=555 + int(th * 1000)) if th > 0 else None
    if eta2 is None:
        rng = np.random.default_rng(7)
        eta2 = (rng.random((T_fft, 8)) > 0.5).astype(np.uint8)
    ml2 = T_fft // 8
    ac_x = np.mean([fft_autocorr(eta2[:, x].astype(np.float64)) for x in range(8)], axis=0)
    tau2 = np.arange(1, ml2)
    s_x = slope_of(tau2, ac_x[1:ml2], 64, 4096)
    th_x = (1.0 + s_x) / 2.0

    noise_rows.append((th, 2 * th - 1, s_f, th_f, p1, s_x, th_x))
    print(f"{th:7.4f} {2*th-1:7.3f} | {s_f:10.3f} {th_f:7.4f} {p1:6.4f} "
          f"| {s_x:10.3f} {th_x:7.4f} | {th_f-th_x:+9.4f}")

print("-" * 104)
print("th_eff = (1 + measured slope)/2, i.e. the exponent a naive asymptotic fit")
print("would report.  A systematic offset between generators at fixed nominal")
print("theta would shift the apparent threshold and is a candidate explanation")
print("for the published spread theta_c = 0 (none) / 1/6 / 0.2 / 0.25.")

# --------------------------------------------------------------------------
# B. theta = 0 dynamics
# --------------------------------------------------------------------------
print()
print("=" * 104)
print("B. theta = 0 DYNAMICS (expect beta = 1/3, alpha = 1/2)")
print("=" * 104)

for rule, rname in [(RULE_BD, "BD"), (RULE_RSOS, "RSOS")]:
    print(f"\n  --- rule = {rname} ---")
    wsat = {}
    for L in [128, 256, 512, 1024]:
        T = 200_000
        spec = SimSpec(theta=0.0, L=L, T=T, n_runs=48, rule=rule, n_check=40, n_w=20000)
        res = simulate(spec)
        wt = res["w_times"].astype(np.float64)
        W = np.sqrt(res["w_series"].astype(np.float64).mean(axis=0))
        wm = W
        wsat[L] = wm[-1]
        b = [slope_of(wt, None, lo, hi) if False else
             np.polyfit(np.log(wt[(wt >= lo) & (wt <= hi)]),
                        np.log(wm[(wt >= lo) & (wt <= hi)]), 1)[0]
             for lo, hi in [(100, 1000), (1000, 10000), (10000, 100000)]]
        print(f"    L={L:5d}  W_sat={wm[-1]:8.3f}  W_sat/sqrt(L)={wm[-1]/np.sqrt(L):.4f}"
              f"   beta[1e2,1e3]={b[0]:.3f} [1e3,1e4]={b[1]:.3f} [1e4,1e5]={b[2]:.3f}"
              f"   ({res['wall_time']:.1f}s)")
    Ls = np.array(sorted(wsat))
    Ws = np.array([wsat[L] for L in Ls])
    al = np.polyfit(np.log(Ls), np.log(Ws), 1)[0]
    print(f"    -> alpha_fit = {al:.4f}   (expect 0.5)")

print("\nDONE")
