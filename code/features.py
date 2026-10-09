"""
Dimensionless, amplitude-free features extracted from interface snapshots.

DESIGN CONSTRAINT (learned the hard way).  For any theta > 0 the noise spectral
density D(omega) = D0 + D_theta omega^{-2 theta} diverges at low frequency, so
the SATURATED WIDTH AMPLITUDE grows with theta even when the universality class
is unchanged.  Any feature carrying that amplitude (e.g. C(1) = <(dh)^2>, the
variance of the height increments) therefore produces a spurious "transition at
theta = 0+" for an anomaly detector.  Every feature below is consequently
normalised to be blind to that amplitude:

  * exponent-like quantities        alpha_loc(r)          (log-log slopes)
  * normalised roughness shape      log10[ C(r) / C(1) ]  (amplitude cancels)
  * normalised distribution shape   standardised moments of dh and of h
  * dimensionless fractions         facet / flat / step fractions, entropy
  * log10(mean flat-run length)     a pure lattice-scale morphology descriptor

The same normalisation is applied to the Hellinger observable, which is computed
on the standardised height-INCREMENT (slope) density rather than on the height
density: the marginal height density turns out to be almost insensitive to
theta (see analyze.py, which reports it as a controlled null result), whereas
the increment density is the quantity that controls faceting.
"""

from __future__ import annotations

import numpy as np

R_LIST = (1, 2, 4, 8, 16, 32, 64, 128, 256)
ALPHA_R = (2, 4, 8, 16, 32, 64, 128, 256)


# --------------------------------------------------------------------------
# basic correlation objects
# --------------------------------------------------------------------------

def height_difference_correlation(h, r_list=R_LIST):
    """C(r) = <(h(x+r) - h(x))^2>, periodic.  h: (..., L) -> (..., len(r_list))."""
    out = np.empty(h.shape[:-1] + (len(r_list),), dtype=np.float64)
    for k, r in enumerate(r_list):
        d = h - np.roll(h, -r, axis=-1)
        out[..., k] = np.mean(d * d, axis=-1)
    return out


def structure_factor(h):
    """S(k) = |FFT(h - mean)|^2 / L for m = 1..L/2 (wavelength L/m sites).

    h : (..., L) -> (..., L//2).
    """
    hh = h - h.mean(axis=-1, keepdims=True)
    F = np.fft.rfft(hh, axis=-1)
    S = (F.real ** 2 + F.imag ** 2) / h.shape[-1]
    return S[..., 1:h.shape[-1] // 2 + 1]


def local_roughness_exponent(C, r_list=R_LIST):
    """alpha_loc(r) = (1/2) d log C / d log r by central differences."""
    lr = np.log(np.asarray(r_list, dtype=np.float64))
    lC = np.log(np.maximum(C, 1e-300))
    al = np.full(C.shape, np.nan, dtype=np.float64)
    for k in range(1, len(r_list) - 1):
        al[..., k] = 0.5 * (lC[..., k + 1] - lC[..., k - 1]) / (lr[k + 1] - lr[k - 1])
    return al


def spectral_roughness_exponent(S, L):
    """alpha_s from S(k) ~ k^{-(2 alpha_s + 1)} fitted over a PHYSICAL window of
    wavelengths 4 to 32 lattice sites, i.e. m/L in [1/32, 1/4].  Using a window
    in m (not in k-fraction) keeps the probed length scale comparable across L.
    """
    n_k = S.shape[-1]
    lo = max(1, int(round(n_k / 16.0)))
    hi = max(lo + 3, int(round(n_k / 4.0)))
    hi = min(hi, n_k)
    m = np.arange(1, n_k + 1, dtype=np.float64)[lo:hi]
    lk = np.log(2.0 * np.pi * m / L)
    lS = np.log(np.maximum(S[..., lo:hi], 1e-300))
    mk = lk.mean()
    denom = ((lk - mk) ** 2).sum()
    slope = ((lS - lS.mean(axis=-1, keepdims=True)) * (lk - mk)).sum(axis=-1) / denom
    return -0.5 * (slope + 1.0)


# --------------------------------------------------------------------------
# shape statistics
# --------------------------------------------------------------------------

def _std_moments(v):
    mu = v.mean(axis=-1, keepdims=True)
    sd = v.std(axis=-1, keepdims=True)
    sd = np.where(sd > 1e-12, sd, 1.0)
    z = (v - mu) / sd
    return z, (z ** 3).mean(axis=-1), (z ** 4).mean(axis=-1) - 3.0


def increment_statistics(h):
    """Statistics of d(x) = h(x+1) - h(x) (periodic), all amplitude-free."""
    d = h - np.roll(h, -1, axis=-1)
    _, sk, ku = _std_moments(d.astype(np.float64))
    b0 = (d == 0)
    starts = b0 & ~np.roll(b0, 1, axis=-1)
    n_runs = np.maximum(starts.sum(axis=-1), 1)
    mean_flat_run = b0.sum(axis=-1) / n_runs

    vals = np.arange(-8, 9)
    flat = d.reshape(-1, d.shape[-1])
    hh = np.empty((flat.shape[0], vals.size), dtype=np.float64)
    for j, v in enumerate(vals):
        hh[:, j] = (flat == v).mean(axis=-1)
    p = hh / np.maximum(hh.sum(axis=-1, keepdims=True), 1e-300)
    ent = (-(p * np.log(np.maximum(p, 1e-300))).sum(axis=-1)
           / np.log(vals.size)).reshape(d.shape[:-1])

    return {
        "skew_d": sk, "kurt_d": ku,
        "frac_flat": b0.mean(axis=-1),
        "frac_step2": (np.abs(d) >= 2).mean(axis=-1),
        "frac_step3": (np.abs(d) >= 3).mean(axis=-1),
        "mean_flat_run": mean_flat_run,
        "slope_entropy": ent,
    }


def height_shape(h):
    hh = h.astype(np.float64)
    _, sk, ku = _std_moments(hh)
    sd = hh.std(axis=-1)
    return {"skew_h": sk, "kurt_h": ku, "W": sd}


# --------------------------------------------------------------------------
# assembled feature vector
# --------------------------------------------------------------------------

def snapshot_features(profiles):
    """Amplitude-free feature matrix for a batch of snapshots.

    profiles : (n_runs, n_check, L) int
    Returns (names, F, extra) with F of shape (n_runs, n_check, n_feat).
    """
    p = profiles.astype(np.float64)
    L = p.shape[-1]
    C = height_difference_correlation(p, R_LIST)
    al = local_roughness_exponent(C, R_LIST)
    S = structure_factor(p)
    a_s = spectral_roughness_exponent(S, L)
    inc = increment_statistics(profiles)
    hs = height_shape(profiles)

    r_idx = {r: i for i, r in enumerate(R_LIST)}
    valid_alpha = [r for r in ALPHA_R if r in r_idx and r < L // 2]
    valid_norm = [r for r in ALPHA_R if r in r_idx and r < L // 2]

    feats = []
    for r in valid_alpha:
        feats.append((f"alpha_loc_r{r}", al[..., r_idx[r]]))
    c1 = np.maximum(C[..., r_idx[1]], 1e-300)
    for r in valid_norm:
        feats.append((f"logCnorm_r{r}", np.log10(np.maximum(C[..., r_idx[r]], 1e-300) / c1)))
    feats += [
        ("alpha_s", a_s),
        ("skew_d", inc["skew_d"]),
        ("kurt_d", inc["kurt_d"]),
        ("frac_flat", inc["frac_flat"]),
        ("frac_step2", inc["frac_step2"]),
        ("frac_step3", inc["frac_step3"]),
        ("slope_entropy", inc["slope_entropy"]),
        ("log_mean_flat_run", np.log10(np.maximum(inc["mean_flat_run"], 1e-6))),
        ("skew_h", hs["skew_h"]),
        ("kurt_h", hs["kurt_h"]),
    ]
    names = [n for n, _ in feats]
    F = np.stack([np.asarray(v, dtype=np.float64) for _, v in feats], axis=-1)
    extra = {"C": C, "alpha_loc": al, "S": S, "alpha_s": a_s, "W": hs["W"],
             "valid_alpha_r": np.asarray(valid_alpha),
             "R_LIST": np.asarray(R_LIST),
             **{k: v for k, v in inc.items()},
             "skew_h": hs["skew_h"], "kurt_h": hs["kurt_h"]}
    return names, F, extra


def alpha_loc_at(extra, r):
    """alpha_loc at a particular separation r (NaN if that r was not computed)."""
    try:
        k = int(np.where(extra["R_LIST"] == r)[0][0])
    except Exception:
        return np.nan
    if k >= extra["alpha_loc"].shape[-1]:
        return np.nan
    return float(np.nanmean(extra["alpha_loc"][..., k]))


# --------------------------------------------------------------------------
# samples for the distributional diagnostic
# --------------------------------------------------------------------------

def pooled_standardized_slopes(profiles_at_t):
    """Pooled standardised height increments d(x) = h(x+1)-h(x).

    Standardising by the ensemble's own mean and variance removes the
    amplitude, so the resulting density is a pure shape observable.
    """
    p = profiles_at_t.astype(np.float64)
    d = (p - np.roll(p, -1, axis=-1)).ravel()
    d = d - d.mean()
    sd = d.std()
    return d / sd if sd > 1e-12 else d


def pooled_standardized_heights(profiles_at_t):
    """Pooled standardised heights h(x) (the insensitive observable, kept for
    the controlled comparison reported in the manuscript)."""
    p = profiles_at_t.astype(np.float64)
    mu = p.mean(axis=-1, keepdims=True)
    sd = p.std(axis=-1, keepdims=True)
    sd = np.where(sd > 1e-12, sd, 1.0)
    return ((p - mu) / sd).ravel()


def width_of(profiles_at_t):
    return profiles_at_t.astype(np.float64).std(axis=-1)
