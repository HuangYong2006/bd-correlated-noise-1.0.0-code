"""
CORRECTED generator calibration: performed on the BINARISED noise.

The model does not use the Gaussian field eta_theta(x,t) directly; it uses the
thresholded variable eta' = 1[eta > 0].  Thresholding is a strong nonlinearity
that compresses correlations,

    C_bin(tau) = (2/pi) arcsin( C_gauss(tau) ),

so the effective theta of a binarised sequence is NOT the effective theta of the
Gaussian sequence it came from.  Since only the binarised field enters the
dynamics, the calibration must be done after thresholding.

Reference family: the binarised EXACT fractional Gaussian noise, whose
correlation is
    C_ref(tau; theta) = (2/pi) arcsin( C_fGn(tau; theta) ),
    C_fGn(tau; theta) = 0.5[(tau+1)^{2H} - 2 tau^{2H} + (tau-1)^{2H}], H=theta+1/2.
This is the unique family that is (i) exactly white noise at theta = 0, since
C_fGn(tau>0; 0) = 0, and (ii) has the correct asymptotic amplitude
2 theta (theta + 1/2).

For each generator we binarise its output, measure C_bin(tau), and fit
theta_eff by weighted least squares against the reference family.

Usage: python code/calibrate_binarised.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
from numba import njit
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from bd_core import ffgn_weights, ffgn_n_components, _u01

T = 1 << 21
NREAL = 6
FIT_LO, FIT_HI = 1, 4096


@njit(cache=True)
def ffgn_seq(T, Wn, rn, seed):
    N = Wn.shape[0]
    X = np.zeros(N)
    out = np.empty(T, dtype=np.float64)
    base = np.uint64(seed + 1) * np.uint64(2654435761)
    for t in range(T):
        z = _u01(base + np.uint64(t)) - 0.5
        s = 0.0
        for n in range(N):
            xn = rn[n] * X[n] + z
            X[n] = xn
            s += Wn[n] * xn
        out[t] = s
    return out


def fft_ac_binary(b):
    x = b.astype(np.float64)
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    return ac / ac[0]


def ref_binarised(tau, theta):
    H = theta + 0.5
    t = tau.astype(np.float64)
    c = 0.5 * ((t + 1) ** (2 * H) - 2 * t ** (2 * H) + np.abs(t - 1) ** (2 * H))
    return (2.0 / np.pi) * np.arcsin(np.clip(c, -1.0, 1.0))


def fit_theta_eff(ac_bin, taus, w):
    def cost(th):
        if th < 0 or th > 0.498:
            return 1e9
        return float(np.sum(w * (ref_binarised(taus, th) - ac_bin[taus]) ** 2))

    r = minimize_scalar(cost, bounds=(0.0, 0.498), method="bounded",
                        options={"xatol": 1e-5})
    return float(r.x)


def ffgn_binary_ac(theta, rng):
    N = ffgn_n_components(T)
    Wn, rn = ffgn_weights(theta, N)
    acs = []
    for i in range(NREAL):
        s = ffgn_seq(T, Wn, rn, seed=int(rng.integers(1, 10 ** 6)))
        acs.append(fft_ac_binary(s > 0.0))
    return np.mean(acs, axis=0)


def spectral_binary_ac(theta, rng):
    freqs = np.fft.rfftfreq(T)
    amp = np.zeros_like(freqs)
    nz = freqs > 0
    amp[nz] = freqs[nz] ** (-theta) if theta != 0.0 else 1.0
    amp[0] = 0.0
    acs = []
    for _ in range(NREAL):
        u = rng.standard_normal(T)
        U = np.fft.rfft(u)
        U *= amp
        s = np.fft.irfft(U, n=T)
        acs.append(fft_ac_binary(s > 0.0))
    return np.mean(acs, axis=0)


def main():
    taus = np.arange(FIT_LO, FIT_HI + 1)
    w = 1.0 / taus
    rng = np.random.default_rng(777)
    nom = [-0.2, -0.1, 0.0, 0.02, 0.05, 0.10, 0.15, 1 / 6, 0.20, 0.25,
           0.30, 0.35, 0.40, 0.45]

    print("=" * 100)
    print("BINARISED calibration: effective theta of the noise actually used by the model")
    print(f"T = 2^{int(np.log2(T))}, {NREAL} realisations, fit over tau in "
          f"[{FIT_LO},{FIT_HI}] with 1/tau weights")
    print("=" * 100)
    print(f"{'theta_nom':>10} {'FFGN th_eff':>13} {'spec th_eff':>13} "
          f"{'FFGN Cb(1)':>11} {'spec Cb(1)':>11} {'ref Cb(1)':>11}")
    print("-" * 100)

    ff_rows, sp_rows = [], []
    for th in nom:
        a_ff = ffgn_binary_ac(th, rng)
        a_sp = spectral_binary_ac(th, rng)
        tf = fit_theta_eff(a_ff, taus, w)
        ts = fit_theta_eff(a_sp, taus, w)
        ff_rows.append((th, tf))
        sp_rows.append((th, ts))
        print(f"{th:10.4f} {tf:13.4f} {ts:13.4f} {a_ff[1]:11.5f} "
              f"{a_sp[1]:11.5f} {ref_binarised(np.array([1]), th)[0]:11.5f}")
    ff_rows = np.array(ff_rows)
    sp_rows = np.array(sp_rows)

    print("-" * 100)
    th_th = 1.0 / 6.0
    print(f"\nRenormalisation-group threshold: theta_th = 1/6 = {th_th:.4f}")
    print(f"  FFGN(theta_nom = 0)     -> effective theta = {ff_rows[2,1]:.4f}"
          f"   {'ABOVE' if ff_rows[2,1] > th_th else 'BELOW'} theta_th")
    print(f"  SPECTRAL(theta_nom = 0) -> effective theta = {sp_rows[2,1]:.4f}"
          f"   {'ABOVE' if sp_rows[2,1] > th_th else 'BELOW'} theta_th")
    print(f"\n  Minimum effective theta reachable by the FFGN family: "
          f"{ff_rows[:,1].min():.4f}")
    print(f"  Nominal theta needed by FFGN to reach effective theta 1/6: "
          f"{'unreachable' if ff_rows[:,1].min() > th_th else 'reachable'}")
    print(f"\n  Offset at theta_nom = 0:  FFGN {ff_rows[2,1]:+.4f}, "
          f"SPECTRAL {sp_rows[2,1]:+.4f}")

    np.savez(os.path.join(ROOT, "data", "generator_calibration_binarised.npz"),
             ffgn=ff_rows, spectral=sp_rows)
    print("\nSaved data/generator_calibration_binarised.npz")


if __name__ == "__main__":
    main()
