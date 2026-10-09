"""
Exact verification that the generator implements the intended power spectrum,
plus a calibration of the NOMINAL theta (the imposed spectral exponent
S(f) ~ f^(-2 theta)) against the EFFECTIVE local decay exponent of the
correlation function measured in the numerically accessible window.

Rationale: the imposed spectrum S(f) = |f|^(-2*theta) has a hard endpoint at the
Nyquist frequency f = 1/2, which produces an O(tau^-1) correction to
C(tau) ~ tau^(2 theta - 1).  Over a finite window the local log-log slope is
therefore systematically steeper than the asymptotic value 2*theta - 1.  Since a
threshold at theta_c = 1/6 is a statement about the asymptotic (continuum)
spectrum, it is essential to know how large this finite-window bias is.

Run:  python code/check_spectrum.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bd_core import _spectral_amplitude

T = 1 << 21
max_lag = 1 << 15
thetas = [0.02, 0.05, 0.10, 0.15, 1 / 6, 0.20, 0.25, 0.30, 0.40, 0.45]


def fft_autocorr(x):
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conjugate(F), n=nfft)[:n]
    return ac / ac[0]


def loglog_slope(tau, c, lo, hi):
    m = (tau >= lo) & (tau <= hi) & (c > 0)
    p = np.polyfit(np.log(tau[m]), np.log(np.abs(c[m])), 1)
    return p[0]


print("=" * 100)
print(f"Exact spectrum check   T = 2^{int(np.log2(T))}, slope window tau in [100, 8000]")
print("=" * 100)
print(f"{'theta':>8} {'2th-1':>8} {'C_th slope':>11} {'C_emp slope':>12} "
      f"{'max|dC|/C':>10} {'th_eff(1+sl)/2':>15}")
print("-" * 100)

rows = []
rng = np.random.default_rng(20240)
for th in thetas:
    amp = _spectral_amplitude(th, T)
    # exact discrete covariance implied by the imposed amplitude
    cth = np.fft.irfft(amp ** 2, n=T)
    cth = cth / cth[0]
    # empirical covariance, averaged over realisations
    nseq = 6
    acl = []
    for _ in range(nseq):
        u = rng.standard_normal(T)
        U = np.fft.rfft(u)
        U *= amp
        s = np.fft.irfft(U, n=T)
        acl.append(fft_autocorr(s))
    ac = np.mean(acl, axis=0)

    tau = np.arange(1, max_lag)
    s_th = loglog_slope(tau, cth[1:max_lag], 100, 8000)
    s_em = loglog_slope(tau, ac[1:max_lag], 100, 8000)
    denom = np.abs(cth[1:max_lag])
    rel = np.max(np.abs(ac[1:max_lag] - cth[1:max_lag])[denom > 1e-4] / denom[denom > 1e-4])
    th_eff = (1 + s_em) / 2.0
    rows.append((th, 2 * th - 1, s_th, s_em, rel, th_eff))
    print(f"{th:8.4f} {2*th-1:8.3f} {s_th:11.3f} {s_em:12.3f} {rel:10.4f} {th_eff:15.4f}")

print("-" * 100)
print("\nInterpretation:")
print("  * 'max|dC|/C' quantifies the agreement between the empirical covariance")
print("    and the covariance implied EXACTLY by the imposed amplitude; it must be")
print("    at the statistical-noise level (~1e-2 or below for this ensemble size).")
print("  * 'C_th slope' vs '2th-1' shows the finite-window bias of the imposed model.")
print("  * 'th_eff' is the theta a naive asymptotic fit would report.")

np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "spectrum_check.npz")
         if os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
         else "spectrum_check.npz",
         rows=np.array(rows))
print("\nDONE")
