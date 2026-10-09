"""
The two independent diagnostics used to locate theta_c.

DIAGNOSTIC 1 -- Hellinger distance (a *distributional* measure).
  It uses ONLY the marginal probability density of the standardised height
  field s = (h - <h>)/W at a fixed observation time.  The distance to the
  theta = 0 (KPZ) reference ensemble is

        H(P,Q) = sqrt( 1 - sum_k sqrt(p_k q_k) )  in [0,1].

  CRUCIAL CALIBRATION: H computed between two *finite* ensembles is strictly
  positive even when P = Q, and its bias grows as the ensemble shrinks or as
  the correlation length grows.  An uncalibrated H therefore produces a
  spurious "onset" that drifts with system size and sample number.  We remove
  it by evaluating the null floor directly: the reference ensemble is split in
  two disjoint halves and the reference-vs-reference distance is measured with
  exactly the same half-ensemble size as the signal.  The reported observable
  is the bias-corrected excess distance

        dH(theta, L) = < H(theta_half, ref_half) - H(ref_half_a, ref_half_b) >

  together with its z-score against the null distribution.

DIAGNOSTIC 2 -- Isolation Forest (a *geometric / multivariate* measure).
  It uses ONLY the dimensionless scaling-and-morphology feature vector of each
  snapshot (local roughness exponents, slope-distribution moments, step and
  facet fractions, entropy), i.e. disjoint information from diagnostic 1,
  which sees only the 1-D height density.  A forest is trained on the
  theta = 0 ensemble, the decision threshold is calibrated at a fixed false
  positive rate alpha on held-out theta = 0 data, and we record

        Phi(theta, L) = fraction of theta snapshots flagged as anomalous,

  again reported as an excess over a calibrated null and as a z-score.

Both observables have the form  X(theta_c) ~ noise  and  X(theta_c) > noise,
so both are used identically: extract the scale-dependent onset theta*(L) and
extrapolate it.  Agreement between two such different statistical machinery is
the central evidence of this work.
"""

from __future__ import annotations

import numpy as np
from sklearn.ensemble import IsolationForest

# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------


def pooled_standardized(profiles_at_t, mode="diff", r=1):
    """Pooled standardised samples of a 1-D marginal of the interface field.

    mode = "diff"   -> height differences over a separation r,
                       v(x) = h(x+r) - h(x), standardised by the ensemble rms.
                       This is the primary distributional observable: it is a
                       rich, continuous density whose shape encodes the
                       correlation structure at the probed scale.
    mode = "incr"   -> nearest-neighbour increments (r = 1).  Reported as a
                       WEAK diagnostic: because the 0/1-noise BD rule is
                       single-step (|Delta h| <= 1), this density has support
                       on {-1,0,+1} only and is a one-parameter family, so the
                       Hellinger distance built on it saturates immediately.
    mode = "height" -> the heights h(x) themselves.  Reported as a BLIND
                       control: it is essentially insensitive to theta.

    Standardising by the ensemble's OWN rms removes the non-universal
    amplitude growth produced by the low-frequency divergence of D(omega), so
    what remains is a pure shape observable.
    """
    p = profiles_at_t.astype(np.float64)
    if mode == "height":
        mu = p.mean(axis=-1, keepdims=True)
        sd = p.std(axis=-1, keepdims=True)
        sd = np.where(sd > 1e-12, sd, 1.0)
        v = ((p - mu) / sd).ravel()
    elif mode in ("incr", "diff"):
        rr = 1 if mode == "incr" else int(r)
        d = p - np.roll(p, -rr, axis=-1)
        # Standardise EACH RUN by its own rms before pooling.  Without this the
        # pooled density is dominated by run-to-run differences in the
        # (non-universal) fluctuation amplitude, and the split-half null floor
        # of the Hellinger distance becomes larger than the physical signal.
        d = d - d.mean(axis=-1, keepdims=True)
        sd = d.std(axis=-1, keepdims=True)
        sd = np.where(sd > 1e-12, sd, 1.0)
        v = (d / sd).ravel()
    else:
        raise ValueError(mode)
    v = v - v.mean()
    s = v.std()
    return v / s if s > 1e-12 else v


def _pdf(samples, edges):
    cnt, _ = np.histogram(samples, bins=edges)
    tot = cnt.sum()
    return cnt.astype(np.float64) / tot if tot > 0 else cnt.astype(np.float64)


def hellinger(p, q):
    """Hellinger distance between two normalised histograms, in [0, 1]."""
    bc = np.sum(np.sqrt(np.maximum(p, 0.0) * np.maximum(q, 0.0)))
    return float(np.sqrt(max(0.0, 1.0 - min(bc, 1.0))))


# --------------------------------------------------------------------------
# EXACT discrete distributional observable (binning-free)
# --------------------------------------------------------------------------
#
# Histogramming a DISCRETE distribution before computing a Hellinger distance
# is a trap.  Standardising the samples rescales the support, so the spikes of
# two ensembles land in *disjoint* bins; the Bhattacharyya coefficient then
# collapses to the overlap of the empty bins only and H jumps towards its
# maximum irrespective of how small the physical difference is.  We therefore
# work with exact probability mass functions over the integer lattice, computed
# from counts, with no binning and no standardisation.

def increment_pmf(profiles_at_t):
    """Exact pmf of the nearest-neighbour increment, support {-1, 0, +1}.

    profiles_at_t : (n_runs, L).  Returns (n_runs, 3) normalised pmfs ordered as
    (P(Delta h = -1), P(0), P(+1)).
    """
    p = profiles_at_t.astype(np.int64)
    d = p - np.roll(p, -1, axis=-1)
    n = d.shape[-1]
    out = np.empty((d.shape[0], 3), dtype=np.float64)
    out[:, 0] = (d == -1).sum(axis=-1) / n
    out[:, 1] = (d == 0).sum(axis=-1) / n
    out[:, 2] = (d == 1).sum(axis=-1) / n
    return out


def terrace_length_pmf(profiles_at_t, max_len=48):
    """Exact pmf of flat-terrace lengths (runs of Delta h = 0), clipped at max_len.

    Returns (n_runs, max_len + 1) normalised pmfs; last bin collects lengths
    >= max_len.
    """
    p = profiles_at_t.astype(np.int64)
    nr, L = p.shape
    out = np.zeros((nr, max_len + 1), dtype=np.float64)
    for i in range(nr):
        d = p[i] - np.roll(p[i], -1)
        b = (d == 0).astype(np.int8)
        # run lengths of consecutive ones on the periodic lattice
        if b.all():
            out[i, max_len] = 1.0
            continue
        idx = np.flatnonzero(b == 0)          # positions where flats end
        if idx.size == 0:
            out[i, max_len] = 1.0
            continue
        lengths = np.diff(np.append(idx, idx[0] + L)) - 1
        lengths = lengths[lengths > 0]
        if lengths.size == 0:
            out[i, max_len] = 1.0
            continue
        counts = np.bincount(np.minimum(lengths, max_len), minlength=max_len + 1)
        out[i] = counts / counts.sum()
    return out


def _pool_pmf(pmf_stack, idx):
    """Average per-run pmfs over a subset of runs (equal weights)."""
    m = pmf_stack[idx].mean(axis=0)
    s = m.sum()
    return m / s if s > 0 else m


def discrete_hellinger_excess(ref_pmf, th_pmf, n_splits=120, seed=0):
    """Null-calibrated Hellinger distance between two ensembles of exact pmfs.

    ref_pmf, th_pmf : (n_runs, K) per-run probability mass functions.
    The null floor is measured by splitting the reference into two disjoint
    halves of exactly the size used for the signal.
    """
    rng = np.random.default_rng(seed)
    n_ref = ref_pmf.shape[0]
    n_th = th_pmf.shape[0]
    n_half = min(n_ref, n_th) // 2
    if n_half < 4:
        return dict(dH=np.nan, z=np.nan)

    hs, hn = [], []
    for _ in range(n_splits):
        perm = rng.permutation(n_ref)
        pa = _pool_pmf(ref_pmf, perm[:n_half])
        pb = _pool_pmf(ref_pmf, perm[n_half:2 * n_half])
        jt = rng.permutation(n_th)[:n_half]
        pt = _pool_pmf(th_pmf, jt)
        hn.append(hellinger(pa, pb))
        hs.append(hellinger(pt, pa))
    hs = np.asarray(hs)
    hn = np.asarray(hn)
    null_mean = float(hn.mean())
    null_std = float(hn.std(ddof=1)) if hn.size > 1 else 0.0
    sig = float(hs.mean())
    d = sig - null_mean
    return dict(dH=float(d), dH_err=float(hs.std(ddof=1)), Hsig=sig,
                Hnull_mean=null_mean, Hnull_std=null_std,
                z=float(d / null_std) if null_std > 0 else np.inf,
                H_full=hellinger(_pool_pmf(ref_pmf, np.arange(n_ref)),
                                 _pool_pmf(th_pmf, np.arange(n_th))))


def _auto_edges(ref, n_bins):
    lo, hi = np.percentile(ref, [0.1, 99.9])
    pad = 0.05 * (hi - lo)
    return np.linspace(lo - pad, hi + pad, n_bins + 1)


def _split_runs(n_runs, rng):
    idx = rng.permutation(n_runs)
    h = n_runs // 2
    return idx[:h], idx[h:2 * h]


# --------------------------------------------------------------------------
# diagnostic 1: bias-corrected Hellinger distance
# --------------------------------------------------------------------------


def hellinger_excess(ref_prof, th_prof, n_bins=64, n_splits=120, seed=0,
                     edges=None, mode="diff", r=1):
    """Bias-corrected Hellinger distance of a standardised marginal density.

    ref_prof, th_prof : (n_runs, L) profiles at the same time and same L.
    mode : "diff" (primary), "incr" (weak), or "height" (blind control).
    Returns dict with dH, z, and the raw quantities.
    """
    rng = np.random.default_rng(seed)
    n_ref = ref_prof.shape[0]
    n_th = th_prof.shape[0]
    n_half = min(n_ref, n_th) // 2
    if n_half < 4:
        return dict(dH=np.nan, z=np.nan)

    ref_full = pooled_standardized(ref_prof, mode, r)
    if edges is None:
        edges = _auto_edges(ref_full, n_bins)

    hs, hn = [], []
    for _ in range(n_splits):
        perm = rng.permutation(n_ref)
        ra = ref_prof[perm[:n_half]]
        rb = ref_prof[perm[n_half:2 * n_half]]
        jt = rng.permutation(n_th)[:n_half]
        ta = th_prof[jt]

        pa = _pdf(pooled_standardized(ra, mode, r), edges)
        pb = _pdf(pooled_standardized(rb, mode, r), edges)
        pt = _pdf(pooled_standardized(ta, mode, r), edges)
        hn.append(hellinger(pa, pb))
        hs.append(hellinger(pt, pa))

    hs = np.asarray(hs)
    hn = np.asarray(hn)
    null_mean = float(hn.mean())
    null_std = float(hn.std(ddof=1)) if hn.size > 1 else 0.0
    sig_mean = float(hs.mean())
    dH = sig_mean - null_mean
    z = (sig_mean - null_mean) / null_std if null_std > 0 else np.inf

    H_full = hellinger(_pdf(ref_full, edges),
                       _pdf(pooled_standardized(th_prof, mode, r), edges))
    return dict(dH=float(dH), dH_err=float(hs.std(ddof=1)),
                Hsig=sig_mean, Hnull_mean=null_mean, Hnull_std=null_std,
                z=float(z), H_full=float(H_full), n_bins=n_bins, mode=mode, r=r)


def hellinger_null_curve(ref_prof, n_bins=64, n_splits=120, seed=1,
                         mode="diff", r=1):
    """Pure null floor of the Hellinger distance at this (L, n_runs, mode)."""
    rng = np.random.default_rng(seed)
    n_ref = ref_prof.shape[0]
    n_half = n_ref // 2
    edges = _auto_edges(pooled_standardized(ref_prof, mode, r), n_bins)
    vals = []
    for _ in range(n_splits):
        perm = rng.permutation(n_ref)
        ra = ref_prof[perm[:n_half]]
        rb = ref_prof[perm[n_half:2 * n_half]]
        vals.append(hellinger(_pdf(pooled_standardized(ra, mode, r), edges),
                              _pdf(pooled_standardized(rb, mode, r), edges)))
    v = np.asarray(vals)
    return float(v.mean()), float(v.std(ddof=1))


# --------------------------------------------------------------------------
# diagnostic 2: Isolation Forest anomaly fraction
# --------------------------------------------------------------------------


class IForestDetector:
    """Isolation-Forest novelty detector calibrated on the theta = 0 ensemble."""

    def __init__(self, n_estimators=400, alpha=0.02, train_frac=0.6,
                 cal_frac=0.2, seed=0, standardize=True):
        self.n_estimators = n_estimators
        self.alpha = alpha
        self.train_frac = train_frac
        self.cal_frac = cal_frac
        self.seed = seed
        self.standardize = standardize
        self.mu = None
        self.sd = None
        self.forest = None
        self.threshold = None
        self.keep = None

    def _prep(self, X):
        X = np.asarray(X, dtype=np.float64)
        if self.standardize:
            X = (X - self.mu) / self.sd
        return X

    def fit(self, X_ref):
        X_ref = np.asarray(X_ref, dtype=np.float64)
        n = X_ref.shape[0]

        # drop degenerate features (essentially zero variance in the reference)
        sd0 = X_ref.std(axis=0)
        self.keep = sd0 > 1e-8 * max(sd0.max(), 1e-12)
        if self.keep.sum() < 3:
            self.keep = np.ones(X_ref.shape[1], dtype=bool)
        X_ref = X_ref[:, self.keep]

        rng = np.random.default_rng(self.seed)
        perm = rng.permutation(n)
        n_tr = int(self.train_frac * n)
        n_ca = int(self.cal_frac * n)
        tr = X_ref[perm[:n_tr]]
        ca = X_ref[perm[n_tr:n_tr + n_ca]]

        self.mu = tr.mean(axis=0)
        self.sd = tr.std(axis=0)
        self.sd = np.where(self.sd > 1e-12, self.sd, 1.0)

        self.forest = IsolationForest(
            n_estimators=self.n_estimators,
            max_samples=min(256, tr.shape[0]),
            contamination="auto",
            random_state=self.seed,
            n_jobs=1,
        )
        self.forest.fit(self._prep(tr))

        scores_ca = self.forest.score_samples(self._prep(ca))
        # lower score_samples == more anomalous; threshold at the alpha quantile
        self.threshold = float(np.quantile(scores_ca, self.alpha))
        self.cal_scores = scores_ca
        self.cal_mean = float(scores_ca.mean())
        self.cal_std = float(scores_ca.std(ddof=1))
        # null distribution of the flagged fraction, by bootstrapping the
        # calibration set
        boot = []
        rng2 = np.random.default_rng(self.seed + 999)
        for _ in range(400):
            s = rng2.choice(scores_ca, size=scores_ca.size, replace=True)
            boot.append(float((s < self.threshold).mean()))
        self.phi_null_mean = float(np.mean(boot))
        self.phi_null_std = float(np.std(boot, ddof=1))
        self.n_features_used = int(self.keep.sum())
        return self

    def score(self, X):
        X = np.asarray(X, dtype=np.float64)[:, self.keep]
        s = self.forest.score_samples(self._prep(X))
        phi = float((s < self.threshold).mean())
        zmean = (float(s.mean()) - self.cal_mean) / max(self.cal_std, 1e-12)
        return dict(phi=phi, mean_score=float(s.mean()), z_mean=float(zmean),
                    scores=s)


def iforest_excess(X_ref, X_theta, n_estimators=400, alpha=0.02, seed=0,
                   n_boot=40):
    """Calibrated anomalous-fraction excess for one (theta, L, t) cell.

    Returns dPhi = Phi(theta) - <Phi_null> clamped at 0, its z-score, and the
    mean-score z-score, averaged over n_boot random train/calibration splits.
    """
    phis, zs, zms = [], [], []
    for b in range(n_boot):
        det = IForestDetector(n_estimators=n_estimators, alpha=alpha,
                              seed=seed + 1000 * b).fit(X_ref)
        r = det.score(X_theta)
        phis.append(r["phi"] - det.phi_null_mean)
        zs.append((r["phi"] - det.phi_null_mean) / max(det.phi_null_std, 1e-12))
        zms.append(r["z_mean"])
    phis = np.asarray(phis)
    zs = np.asarray(zs)
    zms = np.asarray(zms)
    # sklearn's score_samples is oriented so that HIGHER == more normal; the
    # anomaly signal is therefore the negated standardised mean score.
    return dict(dPhi=float(max(0.0, phis.mean())), dPhi_raw=float(phis.mean()),
                dPhi_err=float(phis.std(ddof=1)),
                z=float(zs.mean()), z_err=float(zs.std(ddof=1)),
                z_anom=float(-zms.mean()), z_anom_err=float(zms.std(ddof=1)),
                z_mean_score=float(zms.mean()), z_mean_score_err=float(zms.std(ddof=1)))


# --------------------------------------------------------------------------
# onset extraction and finite-size scaling of the onset
# --------------------------------------------------------------------------


def onset_theta(thetas, values, threshold, increasing=True):
    """First theta (linearly interpolated) at which ``values`` crosses
    ``threshold``.  Returns nan if there is no crossing inside the grid."""
    th = np.asarray(thetas, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    th, v = th[ok], v[ok]
    if th.size < 2:
        return np.nan
    o = np.argsort(th)
    th, v = th[o], v[o]
    if increasing:
        for i in range(th.size - 1):
            if v[i] < threshold <= v[i + 1]:
                dv = v[i + 1] - v[i]
                if abs(dv) < 1e-15:
                    return float(th[i])
                return float(th[i] + (threshold - v[i]) * (th[i + 1] - th[i]) / dv)
    else:
        for i in range(th.size - 1):
            if v[i] > threshold >= v[i + 1]:
                dv = v[i] - v[i + 1]
                if abs(dv) < 1e-15:
                    return float(th[i])
                return float(th[i] + (v[i] - threshold) * (th[i + 1] - th[i]) / dv)
    return np.nan


def onset_with_error(thetas, values, sigma, threshold, increasing=True,
                     n_boot=400, seed=0):
    """Onset theta* with a bootstrap error bar.

    The observable is perturbed by its own per-point uncertainty and the
    crossing is re-extracted, which propagates the measurement noise of the
    onset curve into sigma(theta*).  Points that never cross are counted as
    failures and reported through ``frac_failed``.
    """
    rng = np.random.default_rng(seed)
    th = np.asarray(thetas, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    s = np.asarray(sigma, dtype=np.float64)
    s = np.where(np.isfinite(s) & (s > 0), s, np.nan)
    central = onset_theta(th, v, threshold, increasing)
    vals = []
    for _ in range(n_boot):
        pert = v + rng.standard_normal(v.shape) * np.where(np.isfinite(s), s, 0.0)
        vals.append(onset_theta(th, pert, threshold, increasing))
    vals = np.asarray(vals, dtype=np.float64)
    good = np.isfinite(vals)
    if good.sum() < 10:
        return dict(theta_star=central, theta_star_err=np.nan,
                    frac_failed=1.0 - good.mean())
    return dict(theta_star=central,
                theta_star_err=float(vals[good].std(ddof=1)),
                theta_star_lo=float(np.percentile(vals[good], 16)),
                theta_star_hi=float(np.percentile(vals[good], 84)),
                frac_failed=float(1.0 - good.mean()))


def fit_onset_scaling(Ls, thetas_star, sigma=None, theta_c_free=True):
    """Fit theta*(L) = theta_c + A * L^(-1/nu)   (theta_c_free=True)
       or theta*(L) =            A * L^(-1/nu)   (theta_c free = 0).

    Returns dict(theta_c, theta_c_err, nu, A, chi2, dof, pvalue).
    """
    from scipy.optimize import least_squares
    from scipy.stats import chi2 as chi2dist

    Ls = np.asarray(Ls, dtype=np.float64)
    y = np.asarray(thetas_star, dtype=np.float64)
    if sigma is None:
        sigma = np.full_like(y, max(np.nanstd(y) / 5.0, 1e-3))
    sigma = np.asarray(sigma, dtype=np.float64)
    ok = np.isfinite(y) & np.isfinite(Ls) & (sigma > 0)
    Ls, y, sigma = Ls[ok], y[ok], sigma[ok]
    if Ls.size < 3:
        return dict(theta_c=np.nan, theta_c_err=np.nan, nu=np.nan, A=np.nan,
                    chi2=np.nan, dof=0, pvalue=np.nan, n_points=int(Ls.size))

    def resid(p):
        if theta_c_free:
            tc, nu, A = p
        else:
            nu, A = p
            tc = 0.0
        return (tc + A * Ls ** (-1.0 / nu) - y) / sigma

    best = None
    grid = ([0.0, 1 / 6, 0.2, 0.25] if theta_c_free else [0.0])
    for tc0 in grid:
        for nu0 in (0.5, 1.0, 2.0, 4.0):
            for A0 in (0.02, 0.1, 0.3):
                p0 = [tc0, nu0, A0] if theta_c_free else [nu0, A0]
                try:
                    r = least_squares(resid, p0, method="lm", max_nfev=20000)
                except Exception:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    if best is None:
        return dict(theta_c=np.nan, theta_c_err=np.nan, nu=np.nan, A=np.nan,
                    chi2=np.nan, dof=0, pvalue=np.nan, n_points=int(Ls.size))

    chi2 = 2.0 * best.cost
    dof = max(1, Ls.size - (3 if theta_c_free else 2))
    p = best.x
    if theta_c_free:
        tc, nu, A = p
    else:
        nu, A = p
        tc = 0.0

    # parameter uncertainty from the Jacobian
    try:
        J = best.jac
        cov = np.linalg.inv(J.T @ J)
        perr = np.sqrt(np.abs(np.diag(cov)))
        tc_err = float(perr[0]) if theta_c_free else 0.0
        nu_err = float(perr[1]) if theta_c_free else float(perr[0])
    except np.linalg.LinAlgError:
        tc_err, nu_err = np.nan, np.nan

    return dict(theta_c=float(tc), theta_c_err=tc_err, nu=float(nu),
                nu_err=nu_err, A=float(A), chi2=float(chi2), dof=int(dof),
                pvalue=float(1.0 - chi2dist.cdf(chi2, dof)),
                n_points=int(Ls.size))


def collapse_quality(theta_c, nu, omega, thetas, Ls, D):
    """Normalised collapse quality Q in (-inf, 1] for the scaling ansatz
    D(theta,L) = L^-omega F( (theta-theta_c) L^(1/nu) ).

    Q = 1 - (residual variance about a local smooth of the master curve)
            / (variance of the rescaled observable).
    Points with theta < theta_c are dropped (there D must vanish).
    """
    th = np.asarray(thetas, dtype=np.float64)[:, None]
    LL = np.asarray(Ls, dtype=np.float64)[None, :]
    Dm = np.asarray(D, dtype=np.float64)
    X = (th - theta_c) * LL ** (1.0 / nu)
    Y = Dm * LL ** omega
    m = np.isfinite(X) & np.isfinite(Y) & (th >= theta_c)
    x, y = X[m], Y[m]
    if x.size < 6:
        return -np.inf, x, y
    o = np.argsort(x)
    x, y = x[o], y[o]
    # local smooth over a window of k nearest neighbours in X
    k = max(3, min(7, x.size // 4))
    res = np.empty_like(y)
    for i in range(x.size):
        lo = max(0, i - k // 2)
        hi = min(x.size, lo + k)
        lo = max(0, hi - k)
        res[i] = y[i] - y[lo:hi].mean()
    var = y.var()
    if var <= 0:
        return -np.inf, x, y
    Q = 1.0 - float(np.mean(res ** 2) / var)
    return Q, x, y
