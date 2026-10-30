"""Exact covariance audit of this repository's shared-innovation AR variant.

For X_n(t)=r_n X_n(t-1)+z_t, eta(t)=sum_n W_n X_n(t), and stationary
zero-mean innovations with variance sigma_z^2, the finite-N covariance is

    C(tau) = sigma_z^2 sum_nm W_n W_m r_n^tau / (1-r_n r_m).

Independent innovations across modes retain only n=m terms.  The production
variant shares one centred uniform innovation, so its cross terms must remain.
It is not strictly Gaussian, and the Gaussian arcsine law for thresholded
correlations is not exact for it.  No historical implementation is assumed.

At finite N the ultimate tail is exponential.  A fitted coefficient of a fixed
power law is only a finite-window descriptor, particularly at theta=0, where an
infinite-scale kernel approximation permits a log(tau)/tau correction.  It is
not a measurement or recalibration of an RG exponent.  The helper below is
stationary; ffgn_seq retains the simulator's zero initial condition.

Usage: python code/check_generator_amplitude.py
The audit prints deterministic results and preserves archived amplitude files.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bd_core import ffgn_weights, ffgn_n_components, _u01


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


def ar_sum_covariance(Wn, rn, lags, *, shared=True, innovation_variance=1.0 / 12.0):
    """Stationary covariance at nonnegative integer lags for stable AR modes.

    ``shared=False`` uses independent innovations of the same variance in each
    mode.  Weights and rates are one-dimensional, with abs(r_n)<1.  The result
    is unnormalised; divide by its zero-lag entry when comparing correlations.
    """
    weights = np.asarray(Wn, dtype=np.float64)
    rates = np.asarray(rn, dtype=np.float64)
    if shared:
        mode_amplitudes = (np.outer(weights, weights)
                           / (1.0 - np.outer(rates, rates))).sum(axis=1)
    else:
        mode_amplitudes = weights ** 2 / (1.0 - rates ** 2)
    return innovation_variance * (rates ** np.asarray(lags)[..., None]
                                  @ mode_amplitudes)


def fit_amplitude(ac, theta, lo, hi):
    """Positive finite-window coefficient at a fixed exponent, not an asymptote.

    A signed or unresolved window returns NaN.  Rectifying a Monte Carlo noise
    floor, or selecting only its positive points, would manufacture amplitude.
    """
    tau = np.arange(1, ac.size)
    m = (tau >= lo) & (tau <= hi)
    if m.sum() < 5:
        return np.nan
    expo = 2 * theta - 1
    y = ac[1:][m]
    if not np.all(np.isfinite(y) & (y > 0)):
        return np.nan
    # A = geometric mean of C(tau)/tau^expo
    return float(np.exp(np.mean(np.log(y) - expo * np.log(tau[m]))))


def exact_fgn_amplitude(theta):
    """Asymptotic coefficient of C_fGn(tau)/C_fGn(0) ~ A tau^(2H-2)."""
    H = theta + 0.5
    return H * (2 * H - 1.0)


def spectral_family_amplitude(theta, T=1 << 20, lo=64, hi=4096):
    """Finite-window coefficient of a DC-removed, band-limited power-law spectrum.

    This is not the exact-fGn spectrum or an asymptotic amplitude.  In particular
    at theta=0 the removed DC mode gives negative finite-record off-zero values,
    so there is no positive coefficient to report.
    """
    freqs = np.fft.rfftfreq(T)
    amp = np.zeros_like(freqs)
    nz = freqs > 0
    amp[nz] = freqs[nz] ** (-theta) if theta != 0.0 else 1.0
    c = np.fft.irfft(amp ** 2, n=T)
    c = c / c[0]
    return fit_amplitude(c, theta, lo, hi)


def main():
    T = 1 << 20
    lo, hi = 64, 4096
    lags = np.arange(hi + 1)
    thetas = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.1667, 0.25, 0.30, 0.40, 0.45]

    print("Repository shared-innovation variant: exact stationary finite-N covariance")
    print(f"N={ffgn_n_components(T)}; positive-window descriptor over tau=[{lo},{hi}]")
    print("A_window fixes the nominal exponent; A_fGn is an asymptotic reference only.")
    print(f"{'theta':>7} {'C1_shared':>12} {'C1_indep':>12} "
          f"{'A_window':>12} {'A_fGn':>12}")

    for th in thetas:
        N = ffgn_n_components(T)
        Wn, rn = ffgn_weights(th, N)
        shared = ar_sum_covariance(Wn, rn, lags)
        independent = ar_sum_covariance(Wn, rn, lags, shared=False)
        shared /= shared[0]
        independent /= independent[0]
        window = fit_amplitude(shared, th, lo, hi)
        print(f"{th:7.4f} {shared[1]:12.6f} {independent[1]:12.6f} "
              f"{window:12.6f} {exact_fgn_amplitude(th):12.6f}")

    probe = np.array([64, 256, 1024, 4096])
    print("\ntheta=0 raw AR family: tau*C(tau)/C(0); a pure 1/tau tail would be flat")
    print(f"{'N':>4} " + " ".join(f"{'tau=' + str(t):>12}" for t in probe))
    for N in (16, 23, 32, 48):
        Wn, rn = ffgn_weights(0.0, N)
        cov = ar_sum_covariance(Wn, rn, np.r_[0, probe])
        print(f"{N:4d} " + " ".join(f"{v:12.6f}" for v in probe * cov[1:] / cov[0]))
    print("\nFinite-N convergence in this window is not a proof of an infinite-lag power law.")
    print("The simulator uses a separate i.i.d. shortcut at theta=0; this table evaluates")
    print("the unbranched AR construction. No RG phase or historical-code inference is made.")


if __name__ == "__main__":
    main()
