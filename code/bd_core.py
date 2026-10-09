"""
Ballistic deposition (and RSOS) with long-range temporally correlated noise.

This module implements the *literature-standard* discrete model, so that the
results are directly comparable with the published numbers for theta_c.  The
model is (Ales & Lopez 2019; Song & Xia 2021; Wang & Xia 2025; Jiang & Xia 2026)

    h(x, t+1) = max[ h(x,t) + eta'_theta(x,t),  h(x-1,t),  h(x+1,t) ]      (BD)
    h(x, t+1) = h(x,t) + eta'_theta(x,t)   if |h(x,t)-h(x+-1,t)| <= 1      (RSOS)
                                           else  h(x,t)                    (rejected)

with a PARALLEL update, flat initial condition and periodic boundaries.  The
noise is a BINARISED, temporally long-range correlated sequence,

    eta_theta(x,t)   Gaussian-like,  <eta(x,t) eta(x,t+tau)> ~ tau^(2*theta-1),
                     statistically independent between different sites x;
    eta'_theta(x,t)  =  1  if eta_theta(x,t) > 0  else  0        (P = 1/2),

which is the discrete realisation of the continuum correlator
D(omega,k) = D0 + D_theta * omega^(-2 theta) (spatially delta-correlated,
temporally power-law correlated).  theta = 0 is ordinary uncorrelated BD.

TWO INDEPENDENT NOISE GENERATORS are provided, because the effective
correlation exponent at a given nominal theta depends on the generator:

  1. ``ffgn_weights``  -- Mandelbrot's Fast Fractional Gaussian Noise
     (the generator used by Ales & Lopez, Song & Xia, Jiang & Xia, Wang & Xia).
     Streaming, O(N) state per site, so arbitrarily long runs are cheap in
     memory.  This is the DEFAULT and the literature-comparable choice.
  2. ``fft_noise_field`` -- exact power-law spectral filtering plus Gaussian
     thresholding, included as a generator-robustness cross-check.

References for the model definition are collected in
``LR_temporal_correlation_numerics_report.md`` in the workspace root.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, asdict

import numpy as np
from numba import njit, prange

RULE_BD = 0
RULE_RSOS = 1

# Mandelbrot FFGN constants (as used throughout the literature)
B_FFGN = 2.0
A_FFGN = 6.0


def ffgn_n_components(T: int) -> int:
    """Number of FFGN components resolving correlation times up to ~T.

    Mandelbrot's rule: the largest scale u_1 = a/B = 3 corresponds to the
    shortest resolved correlation time, and u_N = a*B^-N must reach ~1/T,
    i.e. N ~ log_B(a*T).  Four extra components are added as a safety margin.
    """
    return int(max(16, np.ceil(np.log(A_FFGN * T) / np.log(B_FFGN)) + 4))


def ffgn_weights(theta: float, N: int):
    """Weights W_n and AR(1) parameters r_n of the FFGN (Mandelbrot).

        u_n = a B^-n,   r_n = exp(-u_n),
        W_n = sqrt( (1-r_n^2) (B^{1/2-theta} - B^{theta-1/2}) u_n^{1-2theta}
                    / Gamma(2-2 theta) )
    """
    from math import gamma
    n = np.arange(1, N + 1, dtype=np.float64)
    u = A_FFGN * B_FFGN ** (-n)
    r = np.exp(-u)
    shape = B_FFGN ** (0.5 - theta) - B_FFGN ** (theta - 0.5)
    if shape <= 0.0:
        # happens only at theta >= 1/2, where FFGN degenerates
        W = np.zeros(N, dtype=np.float64)
    else:
        W = np.sqrt((1.0 - r ** 2) * shape * u ** (1.0 - 2.0 * theta) / gamma(2.0 - 2.0 * theta))
    return W.astype(np.float64), r.astype(np.float64)


# --------------------------------------------------------------------------
# counter-based uniform RNG (splitmix64): deterministic, thread-safe,
# reproducible independently of the number of threads
# --------------------------------------------------------------------------

@njit(inline="always", cache=True)
def _u01(k):
    z = np.uint64(k) + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z = z ^ (z >> np.uint64(31))
    return np.float64(z >> np.uint64(11)) * (1.0 / 9007199254740992.0)


# --------------------------------------------------------------------------
# core ensemble simulator
# --------------------------------------------------------------------------

@njit(parallel=True, cache=True, fastmath=True, nogil=True)
def _simulate_core(seeds, L, T, rule, theta, Wn, rn, stride, check_idx,
                   profiles, wseries, out_stats):
    """One independent realisation per entry of ``seeds``, parallel over runs.

    profiles : (n_runs, n_check, L) int32   height profiles at check_idx
    wseries  : (n_runs, n_w) float32        interface width at strided times
    out_stats: (n_runs, 3) float64          [mean eta', fraction changed, final W]
    """
    n_runs = seeds.shape[0]
    N = Wn.shape[0]
    n_check = check_idx.shape[0]
    n_w = wseries.shape[1]

    for run in prange(n_runs):
        h = np.zeros(L, dtype=np.int64)
        hnew = np.zeros(L, dtype=np.int64)
        X = np.zeros((L, N), dtype=np.float64)
        eta = np.zeros(L, dtype=np.float64)
        base = np.uint64(seeds[run]) * np.uint64(2654435761) + np.uint64(0x2545F4914F6CDD1D)

        c = 0
        iw = 0
        ns = stride
        sum_eta = 0.0
        n_eta = 0

        for t in range(T):
            # ---- generate the binarised long-range correlated noise ----
            if theta == 0.0:
                for x in range(L):
                    k = base + np.uint64(t) * np.uint64(L) + np.uint64(x)
                    eta[x] = 1.0 if _u01(k) > 0.5 else 0.0
            else:
                for x in range(L):
                    k = base + np.uint64(t) * np.uint64(L) + np.uint64(x)
                    z = _u01(k) - 0.5
                    s = 0.0
                    for n in range(N):
                        xn = rn[n] * X[x, n] + z
                        X[x, n] = xn
                        s += Wn[n] * xn
                    eta[x] = 1.0 if s > 0.0 else 0.0

            # ---- parallel surface update ----
            if rule == 0:
                for x in range(L):
                    xm = x - 1
                    if xm < 0:
                        xm = L - 1
                    xp = x + 1
                    if xp >= L:
                        xp = 0
                    m = h[xm]
                    if h[xp] > m:
                        m = h[xp]
                    v = h[x] + np.int64(eta[x])
                    if v > m:
                        m = v
                    hnew[x] = m
            else:
                for x in range(L):
                    xm = x - 1
                    if xm < 0:
                        xm = L - 1
                    xp = x + 1
                    if xp >= L:
                        xp = 0
                    if abs(h[x] - h[xm]) <= 1 and abs(h[x] - h[xp]) <= 1:
                        hnew[x] = h[x] + np.int64(eta[x])
                    else:
                        hnew[x] = h[x]

            for x in range(L):
                h[x] = hnew[x]
                sum_eta += eta[x]
                n_eta += 1

            # ---- bookkeeping ----
            ns -= 1
            if ns == 0:
                ns = stride
                if iw < n_w:
                    s1 = 0.0
                    s2 = 0.0
                    for x in range(L):
                        fh = np.float64(h[x])
                        s1 += fh
                        s2 += fh * fh
                    mu = s1 / L
                    w2 = s2 / L - mu * mu
                    wseries[run, iw] = w2 if w2 > 0.0 else 0.0
                    iw += 1

            if c < n_check and t + 1 == check_idx[c]:
                for x in range(L):
                    profiles[run, c, x] = np.int32(h[x])
                c += 1

        s1 = 0.0
        s2 = 0.0
        for x in range(L):
            fh = np.float64(h[x])
            s1 += fh
            s2 += fh * fh
        mu = s1 / L
        out_stats[run, 0] = sum_eta / n_eta if n_eta > 0 else 0.0
        out_stats[run, 1] = float(c)
        out_stats[run, 2] = np.sqrt(max(s2 / L - mu * mu, 0.0))

    return 0


@dataclass
class SimSpec:
    theta: float
    L: int
    T: int
    n_runs: int
    rule: int = RULE_BD
    n_check: int = 32
    n_w: int = 20000
    seed: int = 12345
    check_frac_min: float = 2.0e-4
    N_ffgn: int = 0          # 0 -> automatic
    spectrum: str = "powerlaw"
    # "powerlaw": |f|^{-theta}, the family used for the campaign reported in the
    #             paper (bare power law, cut off at Nyquist; theta axis accurate
    #             to ~10%).
    # "fgn":      the exact fractional-Gaussian-noise spectral density, which is
    #             validated in code/validate_fgn_spectrum.py to give a theta axis
    #             accurate to ~1% and to be exactly white at theta = 0.  Both
    #             reduce to white noise at theta = 0.

    def check_times(self):
        t0 = int(max(16, round(self.T * self.check_frac_min)))
        ts = np.unique(np.round(np.geomspace(t0, self.T, self.n_check)).astype(np.int64))
        ts = ts[(ts >= 1) & (ts <= self.T)]
        return ts


def simulate(spec: SimSpec, verbose=False):
    ts = spec.check_times()
    n_check = len(ts)
    stride = max(1, spec.T // spec.n_w)
    n_w = spec.T // stride
    N = spec.N_ffgn if spec.N_ffgn > 0 else ffgn_n_components(spec.T)
    Wn, rn = ffgn_weights(spec.theta, N)

    profiles = np.zeros((spec.n_runs, n_check, spec.L), dtype=np.int32)
    wseries = np.zeros((spec.n_runs, n_w), dtype=np.float32)
    out_stats = np.zeros((spec.n_runs, 3), dtype=np.float64)

    ss = np.random.SeedSequence(spec.seed)
    seeds = ss.generate_state(spec.n_runs, dtype=np.uint32).astype(np.int64)

    t0 = time.time()
    _simulate_core(seeds, spec.L, spec.T, spec.rule, float(spec.theta),
                   Wn, rn, stride, ts, profiles, wseries, out_stats)
    wall = time.time() - t0

    if verbose:
        print(f"    theta={spec.theta:+.4f} L={spec.L} T={spec.T} "
              f"runs={spec.n_runs} rule={spec.rule} N={N}: {wall:.1f}s "
              f"(eta' mean={out_stats[:,0].mean():.4f})", flush=True)

    return {
        "spec": asdict(spec),
        "N_ffgn": N,
        "check_times": ts,
        "w_series": wseries,
        "w_times": np.arange(1, n_w + 1, dtype=np.int64) * stride,
        "profiles": profiles,
        "stats": out_stats,
        "wall_time": wall,
    }


# --------------------------------------------------------------------------
# alternative generator: exact power-law spectral filtering (cross-check)
# --------------------------------------------------------------------------

def fft_noise_field(theta, L, T, seed):
    """Binarised noise field eta'(x,t) in {0,1} from exact spectral filtering.

    Gaussian white noise -> rFFT -> amplitude |f|^(-theta) -> irFFT gives a
    stationary Gaussian field with S(f) ~ |f|^(-2 theta), hence
    C(tau) ~ tau^(2 theta - 1); thresholding at zero binarises it while
    preserving the power-law long-range dependence.
    Returns (T, L) uint8.  Cost O(L*T log T), so use for modest L*T only.
    """
    freqs = np.fft.rfftfreq(T)
    amp = np.zeros_like(freqs)
    nz = freqs > 0
    amp[nz] = freqs[nz] ** (-theta) if theta != 0.0 else 1.0
    rng = np.random.default_rng(seed)
    out = np.empty((T, L), dtype=np.uint8)
    step = max(1, min(L, int(4e6) // max(T, 1)))
    for x0 in range(0, L, step):
        nb = min(step, L - x0)
        u = rng.standard_normal((nb, T))
        U = np.fft.rfft(u, axis=1)
        U *= amp
        s = np.fft.irfft(U, n=T, axis=1)
        out[:, x0:x0 + nb] = (s.T > 0.0).astype(np.uint8)
    return out


@njit(parallel=True, cache=True, fastmath=True, nogil=True)
def _simulate_from_field(eta, rule, L, T, stride, check_idx, profiles, wseries):
    """Same dynamics but driven by a precomputed noise field eta (T, L) uint8."""
    n_runs = eta.shape[0] if eta.ndim == 3 else 1
    n_check = check_idx.shape[0]
    n_w = wseries.shape[1]
    for run in prange(n_runs):
        h = np.zeros(L, dtype=np.int64)
        hnew = np.zeros(L, dtype=np.int64)
        c = 0
        iw = 0
        ns = stride
        for t in range(T):
            for x in range(L):
                xm = x - 1
                if xm < 0:
                    xm = L - 1
                xp = x + 1
                if xp >= L:
                    xp = 0
                e = np.int64(eta[run, t, x])
                if rule == 0:
                    m = h[xm]
                    if h[xp] > m:
                        m = h[xp]
                    v = h[x] + e
                    if v > m:
                        m = v
                    hnew[x] = m
                else:
                    if abs(h[x] - h[xm]) <= 1 and abs(h[x] - h[xp]) <= 1:
                        hnew[x] = h[x] + e
                    else:
                        hnew[x] = h[x]
            for x in range(L):
                h[x] = hnew[x]
            ns -= 1
            if ns == 0:
                ns = stride
                if iw < n_w:
                    s1 = 0.0
                    s2 = 0.0
                    for x in range(L):
                        fh = np.float64(h[x])
                        s1 += fh
                        s2 += fh * fh
                    mu = s1 / L
                    w2 = s2 / L - mu * mu
                    wseries[run, iw] = w2 if w2 > 0.0 else 0.0
                    iw += 1
            if c < n_check and t + 1 == check_idx[c]:
                for x in range(L):
                    profiles[run, c, x] = np.int32(h[x])
                c += 1
    return 0


# --------------------------------------------------------------------------
# I/O helpers
# --------------------------------------------------------------------------

def save_result(path, res):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(
        path,
        spec=json.dumps(res["spec"]),
        N_ffgn=np.int64(res["N_ffgn"]),
        check_times=res["check_times"],
        w_times=res["w_times"],
        w_series=res["w_series"],
        profiles=res["profiles"],
        stats=res["stats"],
        wall_time=np.float64(res["wall_time"]),
    )


def load_result(path):
    z = np.load(path, allow_pickle=False)
    out = {k: z[k] for k in z.files if k not in ("spec", "N_ffgn", "wall_time")}
    out["spec"] = json.loads(str(z["spec"]))
    out["N_ffgn"] = int(z["N_ffgn"])
    out["wall_time"] = float(z["wall_time"])
    return out


def empirical_correlation(seq, max_lag):
    x = np.asarray(seq, dtype=np.float64)
    x = x - x.mean()
    n = x.size
    ac = np.empty(max_lag, dtype=np.float64)
    for tau in range(max_lag):
        ac[tau] = np.dot(x[:n - tau], x[tau:]) / (n - tau)
    return ac / ac[0]


# --------------------------------------------------------------------------
# theta-continuous generator: exact power-law spectral filtering
# --------------------------------------------------------------------------
#
# The FFGN family used throughout the published literature is NOT continuous at
# theta = 0: its normalised correlation amplitude stays O(1) for every theta
# (see code/calibrate_generator.py), so FFGN(theta=0) is not white noise and no
# member of the family reaches the short-range KPZ regime.  The spectral family
#
#     S(f) = |f|^(-2 theta)   (DC mode removed)
#
# IS continuous: at theta = 0 it is exactly white noise, and its correlation
# amplitude tracks the exact fractional-Gaussian-noise amplitude to within
# ~20% over the whole range.  It is therefore the correct generator for a
# theta-scan that is meant to be compared with theta_th = 1/6.
#
# Cost is O(L T log T), so it is used for moderate L; scipy.fft threads and
# 5-smooth FFT lengths keep it affordable.


def smooth_ceil(n: int) -> int:
    """Smallest 5-smooth integer >= n (fast FFT length)."""
    best = None
    p2 = 1
    while p2 < 4 * n:
        p3 = p2
        while p3 < 4 * n:
            p5 = p3
            while p5 < 4 * n:
                if p5 >= n and (best is None or p5 < best):
                    best = p5
                p5 *= 5
            p3 *= 3
        p2 *= 2
    return int(best if best is not None else n)


def _amp(theta, T):
    f = np.fft.rfftfreq(T)
    a = np.zeros_like(f)
    nz = f > 0
    a[nz] = f[nz] ** (-theta) if theta != 0.0 else 1.0
    a[0] = 0.0
    return a


def fgn_spectral_amp(theta, T, n_terms=40):
    """Amplitude of the EXACT fractional-Gaussian-noise spectral density.

    S(f) = (1 - cos 2 pi f) * sum_k |2 pi (f + k)|^{-(2H+1)},   H = theta + 1/2,

    evaluated on the rFFT frequency grid (DC removed).  Near f -> 0 this behaves
    as |f|^{1-2H} = |f|^{-2 theta}, so it realises the required long-range
    exponent, but it also carries the correct high-frequency behaviour and the
    exact correlation AMPLITUDE.  A bare power law |f|^{-2 theta} is cut off at
    the Nyquist frequency and consequently under-correlates at large theta
    (A/A_fGn = 0.61 at theta = 0.4).  This family reduces to white noise at
    theta = 0 because S(f) -> 2(1-cos 2 pi f) * const is flat.
    """
    H = theta + 0.5
    f = np.fft.rfftfreq(T)[1:]                    # drop DC
    S = np.zeros_like(f)
    for k in range(-n_terms, n_terms + 1):
        w = np.abs(2.0 * np.pi * (f + k))
        S += w ** (-(2.0 * H + 1.0))
    S *= (1.0 - np.cos(2.0 * np.pi * f))
    amp = np.sqrt(S)
    out = np.zeros(f.size + 1)
    out[1:] = amp
    return out


def _amp_any(theta, T, spectrum="powerlaw"):
    if spectrum == "fgn":
        return fgn_spectral_amp(theta, T)
    return _amp(theta, T)


def generate_spectral_eta(theta, L, T, n_runs, seed, rows_per_chunk=8,
                          max_bytes=3.0e8, spectrum="powerlaw"):
    """Yield (run0, eta) with eta (nb, T, L) uint8, binarised spectral noise.

    One independent field per run; site rows are processed in small chunks so
    that memory stays bounded and numpy/scipy FFT calls stay cache friendly.
    ``spectrum`` selects the filtering amplitude: "powerlaw" for |f|^{-theta}
    or "fgn" for the exact fractional-Gaussian-noise spectral density.
    """
    from scipy import fft as sfft
    per_run = T * L
    nb = int(max(1, min(n_runs, max_bytes // max(per_run, 1))))
    amp = _amp_any(theta, T, spectrum)
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 991]))
    for r0 in range(0, n_runs, nb):
        n = min(nb, n_runs - r0)
        eta = np.empty((n, T, L), dtype=np.uint8)
        flat = eta.reshape(n * L, T)
        for c0 in range(0, n * L, rows_per_chunk):
            c1 = min(c0 + rows_per_chunk, n * L)
            u = rng.standard_normal((c1 - c0, T))
            U = sfft.rfft(u, axis=1, workers=-1)
            U *= amp
            s = sfft.irfft(U, n=T, axis=1, workers=-1)
            flat[c0:c1] = (s > 0.0)
            del u, U, s
        yield r0, eta
        del eta, flat


def simulate_spectral(spec: SimSpec, verbose=False):
    """Same as ``simulate`` but driven by the theta-continuous spectral generator."""
    ts = spec.check_times()
    n_check = len(ts)
    stride = max(1, spec.T // spec.n_w)
    n_w = spec.T // stride

    profiles = np.zeros((spec.n_runs, n_check, spec.L), dtype=np.int32)
    wseries = np.zeros((spec.n_runs, n_w), dtype=np.float32)

    t0 = time.time()
    for r0, eta in generate_spectral_eta(spec.theta, spec.L, spec.T, spec.n_runs,
                                         spec.seed,
                                         spectrum=getattr(spec, "spectrum", "fgn")):
        n = eta.shape[0]
        _simulate_from_field(eta, spec.rule, spec.L, spec.T, stride, ts,
                             profiles[r0:r0 + n], wseries[r0:r0 + n])
    wall = time.time() - t0
    if verbose:
        print(f"    [spectral] theta={spec.theta:+.4f} L={spec.L} T={spec.T} "
              f"runs={spec.n_runs}: {wall:.1f}s", flush=True)
    return {
        "spec": asdict(spec),
        "N_ffgn": 0,
        "generator": f"spectral-{getattr(spec, 'spectrum', 'fgn')}",
        "check_times": ts,
        "w_series": wseries,
        "w_times": np.arange(1, n_w + 1, dtype=np.int64) * stride,
        "profiles": profiles,
        "stats": np.zeros((spec.n_runs, 3)),
        "wall_time": wall,
    }
