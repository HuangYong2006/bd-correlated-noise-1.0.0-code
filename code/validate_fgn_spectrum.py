"""
Validate the exact-fGn-spectrum generator against the bare power law.

A bare power-law filter |f|^{-2 theta} is cut off at the Nyquist frequency, so
it under-correlates relative to exact fractional Gaussian noise as theta grows
(A/A_fGn = 0.61 at theta = 0.4).  That deficit is the main caveat on our
corrected-generator results, because a weaker long-range coupling implies a
larger crossover length and hence a crossover pushed beyond the accessible L.

Filtering with the exact fGn spectral density
    S(f) = (1 - cos 2 pi f) * sum_k |2 pi (f+k)|^{-(2H+1)},  H = theta + 1/2,
should fix both the exponent and the amplitude.

Usage: python code/validate_fgn_spectrum.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from bd_core import _amp, fgn_spectral_amp

T = 1 << 20
NREAL = 5
FIT_HI = 4096


def fft_ac(x):
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    return ac / ac[0]


def gen(amp, rng):
    u = rng.standard_normal(T)
    U = np.fft.rfft(u)
    U *= amp
    return np.fft.irfft(U, n=T)


def ref_bin(tau, theta):
    H = theta + 0.5
    t = tau.astype(np.float64)
    c = 0.5 * ((t + 1) ** (2 * H) - 2 * t ** (2 * H) + np.abs(t - 1) ** (2 * H))
    return (2.0 / np.pi) * np.arcsin(np.clip(c, -1.0, 1.0))


def fit_theta(acb, taus, w):
    def cost(th):
        if th < 0 or th > 0.498:
            return 1e9
        return float(np.sum(w * (ref_bin(taus, th) - acb[taus]) ** 2))
    r = minimize_scalar(cost, bounds=(0.0, 0.498), method="bounded",
                        options={"xatol": 1e-5})
    return float(r.x)


def amp_ratio(amp, theta, rng):
    """Amplitude of C(tau) ~ A tau^(2 theta - 1), relative to 2 theta(theta+1/2)."""
    acs = [fft_ac(gen(amp, rng)) for _ in range(NREAL)]
    ac = np.mean(acs, axis=0)
    tau = np.arange(1, FIT_HI + 1)
    expo = 2 * theta - 1
    y = np.abs(ac[1:FIT_HI + 1])
    ok = y > 1e-12
    A = float(np.exp(np.mean(np.log(y[ok]) - expo * np.log(tau[ok]))))
    exact = 2 * theta * (theta + 0.5)
    return A, (A / exact if exact > 0 else np.nan)


def main():
    rng = np.random.default_rng(31415)
    thetas = [0.0, 0.02, 0.05, 0.10, 0.15, 1 / 6, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    w = 1.0 / np.arange(1, FIT_HI + 1)
    taus = np.arange(1, FIT_HI + 1)

    print("=" * 104)
    print("Correlation amplitude relative to exact fGn  A_fGn = 2*theta*(theta+1/2)")
    print(f"T = 2^{int(np.log2(T))}, {NREAL} realisations")
    print("=" * 104)
    print(f"{'theta':>7} {'A_fGn':>10} | {'powerlaw A':>11} {'ratio':>8} "
          f"| {'fGn-spec A':>11} {'ratio':>8}")
    print("-" * 104)
    for th in thetas:
        if th == 0.0:
            continue
        ap = _amp(th, T)
        af = fgn_spectral_amp(th, T)
        Ap, rp = amp_ratio(ap, th, rng)
        Af, rf = amp_ratio(af, th, rng)
        print(f"{th:7.4f} {2*th*(th+0.5):10.4f} | {Ap:11.4f} {rp:8.3f} "
              f"| {Af:11.4f} {rf:8.3f}")

    print("\n" + "=" * 104)
    print("Binarised effective theta of both spectral families")
    print("=" * 104)
    print(f"{'theta_nom':>10} {'powerlaw theff':>16} {'fGn-spec theff':>16}")
    print("-" * 104)
    for th in thetas:
        acp = np.mean([fft_ac((gen(_amp(th, T), rng) > 0).astype(float))
                       for _ in range(NREAL)], axis=0)
        acf = np.mean([fft_ac((gen(fgn_spectral_amp(th, T), rng) > 0).astype(float))
                       for _ in range(NREAL)], axis=0)
        tp = fit_theta(acp, taus, w)
        tf = fit_theta(acf, taus, w)
        print(f"{th:10.4f} {tp:16.4f} {tf:16.4f}")

    print("\nInterpretation: the fGn-spectrum column should have ratio ~ 1 across")
    print("the whole range, which removes the amplitude caveat.")
    print("theta = 0 must give exactly white noise for both families.")


if __name__ == "__main__":
    main()
