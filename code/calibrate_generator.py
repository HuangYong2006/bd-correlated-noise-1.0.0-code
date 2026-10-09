"""
CALIBRATION OF THE FFGN theta-AXIS.

The FFGN generator used by every published study of this problem produces a
noise whose NORMALISED correlation function C(tau)/C(0) does not match the
exact fractional-Gaussian-noise shape at the same nominal theta.  Because the
discrete model binarises the noise, eta' = 1[eta > 0], and the binarisation map
C_bin(tau) = (2/pi) arcsin(C(tau)) is monotone, ONLY the normalised shape
C(tau)/C(0) matters for the physics.  Hence the physically meaningful parameter
is not the nominal FFGN theta but the EFFECTIVE theta defined by matching the
normalised correlation function to the exact fGn form

    C_fGn(tau)/C_fGn(0) = (1/2)[(tau+1)^{2H} - 2 tau^{2H} + (tau-1)^{2H}],
    H = theta_eff + 1/2.

This script measures C(tau)/C(0) from long FFGN sequences over a grid of
nominal theta and fits theta_eff by least squares, giving the mapping
theta_eff(theta_nominal).  It then inverts it, so that a threshold reported in
nominal FFGN units (as in the published literature) can be translated into the
effective units in which the renormalisation-group prediction theta_th = 1/6 is
stated.

Usage: python code/calibrate_generator.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
from numba import njit
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bd_core import ffgn_weights, ffgn_n_components, _u01

T_SEQ = 1 << 22
N_REAL = 8
MAX_LAG = 1 << 13
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


def fft_ac(x):
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    return ac / ac[0]


def measured_norm_corr(theta_nom):
    N = ffgn_n_components(T_SEQ)
    Wn, rn = ffgn_weights(theta_nom, N)
    acs = [fft_ac(ffgn_seq(T_SEQ, Wn, rn, seed=1000 + i)) for i in range(N_REAL)]
    return np.mean(acs, axis=0)


def exact_norm_corr(tau, theta_eff):
    H = theta_eff + 0.5
    t = tau.astype(np.float64)
    c = 0.5 * ((t + 1.0) ** (2 * H) - 2.0 * t ** (2 * H) + np.abs(t - 1.0) ** (2 * H))
    return c


def fit_theta_eff(ac, taus, weights):
    def cost(th):
        if th < 0.0 or th > 0.499:
            return 1e9
        model = exact_norm_corr(taus, th)
        return float(np.sum(weights * (model - ac[taus]) ** 2))

    res = minimize_scalar(cost, bounds=(0.0, 0.498), method="bounded",
                          options={"xatol": 1e-5})
    return float(res.x), float(res.fun)


def main():
    taus = np.arange(FIT_LO, FIT_HI + 1)
    # log-spaced weights so that all decades contribute equally
    weights = 1.0 / taus

    nom = [0.0, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.075, 0.10, 0.125, 0.15,
           0.1667, 0.20, 0.225, 0.25, 0.30, 0.35, 0.40, 0.45]

    print("=" * 96)
    print("FFGN theta-axis calibration: nominal theta  ->  effective theta")
    print(f"T = 2^{int(np.log2(T_SEQ))} sweeps, {N_REAL} realisations, "
          f"fit over tau in [{FIT_LO}, {FIT_HI}] with 1/tau weights")
    print("=" * 96)
    print(f"{'theta_nom':>10} {'theta_eff':>10} {'C(1)/C(0)':>11} "
          f"{'C(16)/C(0)':>11} {'C(256)/C(0)':>12} {'resid':>10}")
    print("-" * 96)

    pairs = []
    for th in nom:
        ac = measured_norm_corr(th)
        te, r = fit_theta_eff(ac, taus, weights)
        pairs.append((th, te))
        print(f"{th:10.4f} {te:10.4f} {ac[1]:11.5f} {ac[16]:11.5f} "
              f"{ac[256]:12.6f} {r:10.3e}")
    pairs = np.array(pairs)

    print("-" * 96)
    print("\nInverse mapping: nominal theta required to realise a given effective theta")
    print("-" * 96)
    nom_v, eff_v = pairs[:, 0], pairs[:, 1]
    o = np.argsort(eff_v)
    for target in (1.0 / 6.0, 0.20, 0.25):
        if eff_v.min() <= target <= eff_v.max():
            nom_needed = float(np.interp(target, eff_v[o], nom_v[o]))
        else:
            nom_needed = np.nan
        print(f"  theta_eff = {target:.4f}   <->  theta_nominal ~ {nom_needed:.4f}")

    print("\nKey numbers for the controversy:")
    for tn in (0.0, 0.1667, 0.20, 0.25):
        te = float(np.interp(tn, nom_v, eff_v))
        print(f"  nominal theta = {tn:.4f}  corresponds to effective theta = {te:.4f}")
    te_at_quarter = float(np.interp(0.25, nom_v, eff_v))
    print(f"\n  => A threshold reported at nominal theta_c = 0.25 (the value quoted by")
    print(f"     Ales & Lopez 2019, Song & Xia 2021 and Wang & Xia 2025) corresponds to")
    print(f"     an EFFECTIVE theta of about {te_at_quarter:.3f} in the units in which the")
    print(f"     renormalisation-group prediction theta_th = 1/6 is stated.")

    np.savez(os.path.join(HERE, "..", "data", "generator_calibration.npz"),
             pairs=pairs)
    print("\nDONE")


if __name__ == "__main__":
    main()
