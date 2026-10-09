"""
QUANTITATIVE TEST OF THE FFGN GENERATOR'S CORRELATION AMPLITUDE.

Every published study of BD/KPZ growth with temporally correlated noise
(Alés & López 2019; Song & Xia 2021; Hu, Hao & Xia 2023; Wang & Xia 2025;
Jiang & Xia 2026) generates its noise with Mandelbrot's fast fractional
Gaussian noise (FFGN),

    eta(t) = sum_{n=1..N} W_n X_t(u_n),   u_n = a B^-n,  r_n = exp(-u_n),
    X_t = r_n X_{t-1} + (zeta_t - 1/2),
    W_n^2 = [12(1-r_n^2)/Gamma(2-2theta)] (B^{1/2-theta} - B^{theta-1/2}) u_n^{1-2theta},

and then binarises at zero.  The FFGN is a sum of AR(1) components, hence a
LINEAR filter of one white sequence zeta_t with kernel K(u) = sum_n W_n r_n^u;
its correlation function is therefore

    C_FFGN(tau) ~ [sum_n c u_n^{2-2H} e^{-u_n tau}]  ~  A_FFGN * tau^{2H-2},

with  H = theta + 1/2  and a coefficient A_FFGN that is essentially theta
INDEPENDENT (the geometric u-grid makes the integral converge to a constant).

The exact fractional Gaussian noise correlator, by contrast, is

    C_fGn(tau) = (1/2)[(tau+1)^{2H} - 2 tau^{2H} + (tau-1)^{2H}]
               ~ H(2H-1) tau^{2H-2} = theta(theta+1/2) * 2 * tau^{2H-2},

whose amplitude VANISHES as theta -> 0.  So the FFGN over-correlates by a
factor ~ A_FFGN / (2 theta (theta+1/2)) at small theta, and FFGN(theta -> 0+)
is NOT white noise.

This script measures, for each theta:
  * A_FFGN  : amplitude fitted from a long simulated FFGN sequence, C(0)=1
  * A_fGn   : the exact asymptotic amplitude, (theta+1/2)(2 theta)
  * A_spec  : amplitude of the exact power-law spectral family S(f)=|f|^{-2 theta}
and reports the ratio.  It also checks the N-dependence, to show the effect is
not a truncation artefact.

Usage: python code/check_generator_amplitude.py
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


def fit_amplitude(ac, theta, lo, hi):
    """Fit |C(tau)| = A tau^(2 theta - 1) over an intermediate window."""
    tau = np.arange(1, ac.size)
    m = (tau >= lo) & (tau <= hi)
    if m.sum() < 5:
        return np.nan
    expo = 2 * theta - 1
    y = np.abs(ac[1:][m])
    ok = y > 0
    if ok.sum() < 5:
        return np.nan
    # A = geometric mean of C(tau)/tau^expo
    return float(np.exp(np.mean(np.log(y[ok]) - expo * np.log(tau[m][ok]))))


def exact_fgn_amplitude(theta):
    """Asymptotic coefficient of C_fGn(tau)/C_fGn(0) ~ A tau^(2H-2)."""
    H = theta + 0.5
    return H * (2 * H - 1.0)


def spectral_family_amplitude(theta, T=1 << 20, lo=64, hi=4096):
    """Amplitude of the exact power-law spectral family S(f)=|f|^{-2 theta},
    with the DC mode removed, normalised to C(0)=1."""
    freqs = np.fft.rfftfreq(T)
    amp = np.zeros_like(freqs)
    nz = freqs > 0
    amp[nz] = freqs[nz] ** (-theta) if theta != 0.0 else 1.0
    c = np.fft.irfft(amp ** 2, n=T)
    c = c / c[0]
    return fit_amplitude(c, theta, lo, hi)


def main():
    T = 1 << 20
    n_sites = 8
    lo, hi = 64, 4096
    thetas = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.1667, 0.25, 0.30, 0.40, 0.45]

    print("=" * 104)
    print("FFGN correlation amplitude vs the exact fractional-Gaussian-noise amplitude")
    print(f"T = 2^{int(np.log2(T))}, {n_sites} realisations averaged, "
          f"fit window tau in [{lo}, {hi}], C(0)=1 in every case")
    print("=" * 104)
    print(f"{'theta':>7} {'A_FFGN':>10} {'A_fGn':>10} {'A_spec':>10} "
          f"{'A_FFGN/A_fGn':>13} {'A_FFGN/A_spec':>14}")
    print("-" * 104)

    rows = []
    for th in thetas:
        N = ffgn_n_components(T)
        Wn, rn = ffgn_weights(th, N)
        acs = [fft_ac(ffgn_seq(T, Wn, rn, seed=100 + i)) for i in range(n_sites)]
        ac = np.mean(acs, axis=0)
        A_f = fit_amplitude(ac, th, lo, hi)
        A_e = exact_fgn_amplitude(th)
        A_s = spectral_family_amplitude(th, 1 << 20, lo, hi)
        r1 = A_f / A_e if A_e > 0 else np.inf
        r2 = A_f / A_s if A_s and np.isfinite(A_s) and A_s > 0 else np.inf
        rows.append((th, A_f, A_e, A_s, r1, r2))
        print(f"{th:7.4f} {A_f:10.4f} {A_e:10.5f} {A_s:10.4f} {r1:13.1f} {r2:14.2f}")

    print("-" * 104)
    print("\nInterpretation: A_fGn is the exact asymptotic coefficient, and it VANISHES")
    print("as theta -> 0 (white noise).  A_FFGN stays O(1).  The ratio therefore")
    print("diverges at small theta: the FFGN family does not connect continuously to")
    print("uncorrelated noise, and every theta > 0 is already strongly correlated.")

    # ---- N-dependence: show this is not a truncation artefact ----
    print("\n" + "=" * 104)
    print("N-dependence of A_FFGN  (if the effect were a truncation artefact it")
    print("would disappear as N grows and the slowest mode 1/u_N moves below T)")
    print("=" * 104)
    print(f"{'theta':>7} " + " ".join(f"{'N=' + str(N):>10}" for N in (16, 23, 32, 48)))
    print("-" * 104)
    for th in (0.0, 0.01, 0.05, 0.1667, 0.25, 0.40):
        vals = []
        for N in (16, 23, 32, 48):
            Wn, rn = ffgn_weights(th, N)
            ac = np.mean([fft_ac(ffgn_seq(T, Wn, rn, seed=200 + i))
                          for i in range(4)], axis=0)
            vals.append(fit_amplitude(ac, th, lo, hi))
        print(f"{th:7.4f} " + " ".join(f"{v:10.4f}" for v in vals))

    np.savez(os.path.join(HERE, "..", "data", "generator_amplitude.npz"),
             rows=np.array(rows))
    print("\nDONE")


if __name__ == "__main__":
    main()
