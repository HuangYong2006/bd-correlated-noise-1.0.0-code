"""
CRITICAL VALIDATION: does the FFGN generator reproduce true white noise as
theta -> 0, or does it have a finite-correlation discontinuity at theta = 0?

Why this matters.  The production reference ensemble uses EXACT i.i.d. Bernoulli
noise at theta = 0, whereas every theta > 0 uses Mandelbrot's FFGN (the
generator used throughout the published literature).  FFGN is a finite-N
approximation: it represents C(tau) as a sum of N AR(1) components whose
slowest correlation time is a B^-N, which for N ~ 23 exceeds T.  If the weights
do not suppress those slow modes strongly enough, FFGN(theta -> 0+) is NOT white
noise, and the whole "correlated" family is separated from the theta = 0
reference by a generator artefact instead of by physics.

We therefore compare, at identical L, T and ensemble size:
  * exact i.i.d. Bernoulli(1/2)                       (true uncorrelated BD)
  * FFGN with theta = 1e-9, 1e-6, 1e-4, ... 0.01      (the generator family)
and report the single structural order parameter of the 0/1-noise BD rule,
p_flat = frac{Delta h = 0}, together with W_sat and the noise correlation.

Usage: python code/check_generator_limit.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bd_core import SimSpec, simulate, ffgn_weights, ffgn_n_components, _u01
from numba import njit


@njit(cache=True)
def _ffgn_binary(n_sites, T, Wn, rn, seed):
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


def pflat_and_W(theta, L, T, n_runs, seed=4242):
    spec = SimSpec(theta=theta, L=L, T=T, n_runs=n_runs, n_check=2, n_w=2000, seed=seed)
    res = simulate(spec, verbose=False)
    h = res["profiles"][:, -1, :].astype(np.int64)
    d = h - np.roll(h, -1, axis=-1)
    p_flat = float((d == 0).mean())
    p_up = float((d == +1).mean())
    W = float(np.sqrt(res["w_series"].astype(np.float64).mean(axis=0))[-1])
    mx = int(np.abs(d).max())
    return p_flat, p_up, W, mx


def noise_corr(theta, T=1 << 16, n_sites=16, lags=(1, 2, 4, 8, 16, 64, 256, 1024)):
    N = ffgn_n_components(T)
    Wn, rn = ffgn_weights(theta, N)
    if theta == 0.0:
        rng = np.random.default_rng(11)
        e = (rng.random((T, n_sites)) > 0.5).astype(np.float64)
    else:
        e = _ffgn_binary(n_sites, T, Wn, rn, seed=11).astype(np.float64)
    out = []
    for lag in lags:
        c = np.mean([np.corrcoef(e[:-lag, j], e[lag:, j])[0, 1] for j in range(n_sites)])
        out.append(c)
    return np.array(out)


def main():
    L, T, NR = 256, 61440, 128
    print("=" * 92)
    print("Generator-limit test:  is FFGN(theta -> 0+) the same process as i.i.d. noise?")
    print(f"L = {L}, T = {T}, n_runs = {NR}, N_ffgn = {ffgn_n_components(T)}")
    print("=" * 92)
    print(f"{'theta':>10} {'p_flat':>8} {'p_up':>8} {'W_sat':>8} {'max|dh|':>8}  note")
    print("-" * 92)

    rows = []
    for theta in [0.0, 1e-9, 1e-6, 1e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 5e-2]:
        pf, pu, W, mx = pflat_and_W(theta, L, T, NR)
        tag = "i.i.d. shortcut" if theta == 0.0 else ("FFGN" if theta > 1e-8 else "FFGN, theta~0")
        print(f"{theta:10.2e} {pf:8.4f} {pu:8.4f} {W:8.3f} {mx:8d}  {tag}")
        rows.append((theta, pf, pu, W))

    print("-" * 92)
    p_iid = rows[0][1]
    p_ffgn0 = rows[1][1]
    print(f"\np_flat(i.i.d.)            = {p_iid:.4f}")
    print(f"p_flat(FFGN, theta=1e-9)  = {p_ffgn0:.4f}")
    print(f"absolute discontinuity    = {abs(p_ffgn0 - p_iid):.4f}")
    if abs(p_ffgn0 - p_iid) < 0.01:
        print("=> The FFGN family connects continuously to true white noise.")
        print("   The p_flat jump at small theta is therefore PHYSICS, not a")
        print("   generator artefact, and the theta = 0 i.i.d. reference is valid.")
    else:
        print("=> DISCONTINUITY: FFGN(theta->0+) is NOT white noise.  The theta = 0")
        print("   reference ensemble must then be generated with the SAME generator")
        print("   (i.e. FFGN) for the comparison to be meaningful.")

    print("\n" + "=" * 92)
    print("Noise autocorrelation of the binarised sequence (should vanish for theta->0)")
    print("=" * 92)
    lags = (1, 2, 4, 8, 16, 64, 256, 1024)
    print(f"{'theta':>10} " + " ".join(f"{l:>8d}" for l in lags))
    print("-" * 92)
    for theta in [0.0, 1e-9, 1e-4, 1e-3, 2.5e-3, 1e-2]:
        c = noise_corr(theta)
        print(f"{theta:10.2e} " + " ".join(f"{v:8.4f}" for v in c))

    np.savez(os.path.join(HERE, "..", "data", "generator_limit.npz"),
             rows=np.array([(r[0], r[1], r[2], r[3]) for r in rows]))


if __name__ == "__main__":
    main()
