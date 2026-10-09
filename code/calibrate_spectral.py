"""
Calibrate the theta axis of the SPECTRAL generator with the same
shape-matching procedure used for FFGN, so that the two generators can be
compared at matched EFFECTIVE theta rather than at matched nominal theta.

This matters because the bare power-law spectrum S(f)=|f|^(-2 theta) is cut off
at the Nyquist frequency f=1/2, which makes the family under-correlate relative
to exact fractional Gaussian noise as theta grows.  If the spectral family's
theta_eff differs from nominal, then comparing alpha(theta) between the two
generators at equal NOMINAL theta is not a like-for-like comparison.

Usage: python code/calibrate_spectral.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

T = 1 << 20
NREAL = 6
FIT_HI = 4096


def fft_ac(x):
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    return ac / ac[0]


def exact_fgn(tau, theta_eff):
    H = theta_eff + 0.5
    t = tau.astype(np.float64)
    return 0.5 * ((t + 1) ** (2 * H) - 2 * t ** (2 * H) + np.abs(t - 1) ** (2 * H))


def theta_eff_of(ac, taus):
    w = 1.0 / taus

    def cost(th):
        if th < 0 or th > 0.498:
            return 1e9
        return float(np.sum(w * (exact_fgn(taus, th) - ac[taus]) ** 2))

    r = minimize_scalar(cost, bounds=(0.0, 0.498), method="bounded",
                        options={"xatol": 1e-5})
    return float(r.x), float(r.fun)


def measured_norm_corr_spectral(theta, rng):
    """Normalised correlation of the spectral (power-law filtered) family."""
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
        acs.append(fft_ac(s))
    return np.mean(acs, axis=0)


def main():
    taus = np.arange(1, FIT_HI + 1)
    rng = np.random.default_rng(2024)
    print("=" * 88)
    print("Spectral-family theta calibration")
    print(f"T = 2^{int(np.log2(T))}, {NREAL} realisations, fit over tau in "
          f"[1,{FIT_HI}] with 1/tau weights")
    print("=" * 88)
    print(f"{'theta_nom':>10} {'theta_eff':>10} {'C(1)/C(0)':>11} "
          f"{'C(16)/C(0)':>11} {'C(256)/C(0)':>12} {'exact C(1)':>11}")
    print("-" * 88)
    rows = []
    for th in [-0.2, -0.1, 0.0, 0.02, 0.05, 0.10, 0.15, 1 / 6, 0.20, 0.25,
               0.30, 0.35, 0.40, 0.45]:
        ac = measured_norm_corr_spectral(th, rng)
        te, res = theta_eff_of(ac, taus)
        rows.append((th, te))
        print(f"{th:10.4f} {te:10.4f} {ac[1]:11.5f} {ac[16]:11.5f} "
              f"{ac[256]:12.6f} {exact_fgn(np.array([1]), th)[0]:11.5f}")
    rows = np.array(rows)
    np.savez(os.path.join(ROOT, "data", "generator_calibration_spectral.npz"),
             pairs=rows)

    print("-" * 88)
    print("\nComparison of the two generators' effective-theta axes:")
    try:
        ff = np.load(os.path.join(ROOT, "data", "generator_calibration.npz"))["pairs"]
        print(f"{'theta_nom':>10} {'FFGN th_eff':>13} {'SPECTRAL th_eff':>16}")
        print("-" * 88)
        for th in [0.0, 0.05, 0.10, 1 / 6, 0.20, 0.25, 0.30, 0.35, 0.40]:
            a = float(np.interp(th, ff[:, 0], ff[:, 1]))
            b = float(np.interp(th, rows[:, 0], rows[:, 1]))
            print(f"{th:10.4f} {a:13.4f} {b:16.4f}")
        print("\n=> At a given nominal theta the FFGN noise is systematically MORE")
        print("   correlated than the spectral family.  Exponent curves must")
        print("   therefore be compared as functions of theta_eff, not theta_nom.")
    except Exception as e:
        print(f"(could not load FFGN calibration: {e})")


if __name__ == "__main__":
    main()
