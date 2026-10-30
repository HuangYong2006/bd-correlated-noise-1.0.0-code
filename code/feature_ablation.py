"""
Feature ablation of the Isolation-Forest diagnostic (diagnostic 2, zPhi)
against the Hellinger diagnostic (diagnostic 1, zH) for the fgnspec campaign.

WHY THIS SCRIPT EXISTS
  A referee objected that the manuscript describes the two diagnostics as
  using "disjoint information" without evidence, and asked for (a) the
  correlation between the two diagnostics' outputs and (b) an ablation of the
  Isolation Forest: rerun it with the increment/terrace features removed and
  again with the scaling/roughness-geometry features removed, and check
  whether the conclusions change.

WHAT IT DOES
  * reads the archived Hellinger z-scores from
        data/analysis/<tag>_bd/diagnostics.csv     (column "zH")
  * recomputes the Isolation-Forest z-score from the raw snapshots for THREE
    feature sets that are identical in every other respect:
        FULL       every feature returned by features.snapshot_features
        NO_INC     all features except the increment/terrace family
                   (skew_d, kurt_d, frac_flat, frac_step2, frac_step3,
                    slope_entropy, log_mean_flat_run)
        NO_SCALING all features except the scaling / roughness-geometry family
                   (alpha_loc_r*, logCnorm_r*, alpha_s)
    The height-distribution shape features (skew_h, kurt_h) are kept in all
    three variants: they belong to neither of the two ablated families.
  * writes the per-cell table to data/analysis/ablation_<tag>.csv and the
    summary numbers to data/analysis/ablation_<tag>.json
  * prints Pearson/Spearman correlations, the FULL-vs-ablated changes, the
    fixed-theta slopes of zPhi against ln L, and a verdict on whether the
    qualitative scale-dependence pattern survives the ablation.

HOW IT MATCHES code/analyze.py
  Every choice below is copied from analyze.py; the ONLY difference is which
  columns of the feature matrix are retained:
    * cells are discovered with the same filename regex and the same sign
      convention (tok p0d1000 -> +0.1000, m0d0500 -> -0.0500);
    * the reference ensemble is the theta = 0 cell at the SAME L;
    * the feature matrix retains the last ``--checkpoints`` checkpoints per run:
          k0 = max(0, n_check - checkpoints)
          F  = snapshot_features(profiles[:, k0:, :])   -> (n_runs, n_ck, n_f)
      iforest_excess uses X = F.mean(axis=1), one vector per independent run;
    * realisation matching: with --if-runs R (analyze.py's "--if-runs"),
          m = min(R, n_runs_theta, n_runs_ref); both ensembles retain m runs;
    * the per-cell seed is analyze.py's formula
          seed = args.seed + int(round(theta * 1e4)) + 10 ** 6
      and the SAME seed, n_estimators and n_boot are used for all three
      variants, so the comparison is paired and fair;
    * iforest_excess is called unchanged from code/diagnostics.py, so the
      detector construction (60% train / 20% threshold calibration / remaining
      independent null runs, requested alpha = 0.02, per-resplit seeds
      seed + 1000*b) is identical. Actual FPR is measured on the null runs;
    * the auxiliary fixed-theta ln L slopes are descriptive fits,
      applied to zPhi (and to zH) and carry a least-squares and a Monte-Carlo
      standard error.

DEVIATIONS FROM THE ARCHIVED CANONICAL RUN (stated, not hidden)
  * The exact command line that produced data/analysis/fgnspec_bd/diagnostics.csv
    is not recorded anywhere in the workspace (the console echo in
    data/analysis_fgnspec.txt is truncated).  The protocol constants used here
    are analyze.py's defaults plus --checkpoints 4 and --if-runs 96, which is
    the configuration the README documents for the analogous "spec" campaign.
  * Because the archived Isolation-Forest parameters are unknown, this script
    also reports the correlation computed directly against the ARCHIVED zPhi
    column of diagnostics.csv as a cross-check; that number needs no
    recomputation and is labelled "ARCHIVED" everywhere it appears.

Usage:
  python code/feature_ablation.py
  python code/feature_ablation.py --if-boot 8 --if-runs 96
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from bd_core import load_result                       # noqa: E402
from features import snapshot_features                # noqa: E402
from diagnostics import iforest_excess                # noqa: E402

VARIANTS = ("FULL", "NO_INC", "NO_SCALING")

# feature families, exactly as documented in features.snapshot_features
INCREMENT_FEATURES = frozenset((
    "skew_d", "kurt_d", "frac_flat", "frac_step2", "frac_step3",
    "slope_entropy", "log_mean_flat_run",
))

# theta values whose fixed-theta ln L slope is reported explicitly
REPORT_THETAS = (0.10, 0.1667, 0.20, 0.25, 0.30)

# a slope is called significant when it is positive and > 2 least-squares sigma
SIGMA_SIG = 2.0


# --------------------------------------------------------------------------
# feature-family bookkeeping
# --------------------------------------------------------------------------

def is_scaling(name):
    """scaling / roughness geometry family (documented in features.py)."""
    return (name.startswith("alpha_loc_r") or name.startswith("logCnorm_r")
            or name == "alpha_s")


def variant_columns(names, variant):
    """Indices of the columns retained by ``variant`` (order preserved)."""
    keep = []
    for i, n in enumerate(names):
        if variant == "FULL":
            ok = True
        elif variant == "NO_INC":
            ok = n not in INCREMENT_FEATURES
        elif variant == "NO_SCALING":
            ok = not is_scaling(n)
        else:
            raise ValueError(variant)
        if ok:
            keep.append(i)
    return np.asarray(keep, dtype=int)


# --------------------------------------------------------------------------
# cell discovery (copy of analyze.discover, same regex and sign convention)
# --------------------------------------------------------------------------

def discover(tag, rule):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$",
                     os.path.basename(f))
        if not m:
            continue
        tok = m.group("tok")
        sign = -1.0 if tok[0] == "m" else 1.0
        out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = f
    return out


def read_diagnostics(tag, rule):
    """Archived Hellinger results, keyed by (theta, L)."""
    p = os.path.join(ROOT, "data", "analysis", f"{tag}_{rule}", "diagnostics.csv")
    if not os.path.exists(p):
        raise SystemExit(f"missing {p}; run code/analyze.py first")
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        key = (round(float(r["theta"]), 4), int(float(r["L"])))
        out[key] = dict(zH=float(r["zH"]), dH=float(r["dH"]),
                        zPhi_archived=float(r["zPhi"]),
                        zPhi_archived_err=float(r["zPhi_err"]),
                        dPhi_archived=float(r["dPhi"]))
    return out


# --------------------------------------------------------------------------
# statistics helpers
# --------------------------------------------------------------------------

def corr_pair(x, y):
    """Pearson and Spearman correlation over the finite pairs of (x, y)."""
    from scipy.stats import pearsonr, spearmanr
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return dict(n=n, pearson=np.nan, pearson_p=np.nan,
                    spearman=np.nan, spearman_p=np.nan)
    pr = pearsonr(x[m], y[m])
    sr = spearmanr(x[m], y[m])
    return dict(n=n, pearson=float(pr[0]), pearson_p=float(pr[1]),
                spearman=float(sr[0]), spearman_p=float(sr[1]))


def slope_vs_lnL(L, y, sigma=None, n_mc=2000, seed=0):
    """Least-squares slope of y against ln L, with two standard errors.

    se_ls    analytic OLS standard error of the slope.
    se_boot  Monte-Carlo standard error obtained by perturbing every point by
             its own per-cell uncertainty ``sigma`` and refitting; when
             ``sigma`` is None the OLS error is reused so the two agree.
    """
    L = np.asarray(L, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    m = np.isfinite(y) & np.isfinite(L) & (L > 0)
    L, y = L[m], y[m]
    if L.size < 3:
        return dict(n=int(L.size), slope=np.nan, se_ls=np.nan, se_boot=np.nan,
                    t=np.nan, significant=False)
    x = np.log(L)
    xm = x.mean()
    sxx = float(((x - xm) ** 2).sum())
    b = float(((x - xm) * y).sum() / sxx)
    a = float(y.mean() - b * xm)
    resid = y - (a + b * x)
    dof = L.size - 2
    s2 = float((resid ** 2).sum() / dof) if dof >= 1 else np.nan
    se_ls = float(np.sqrt(s2 / sxx)) if np.isfinite(s2) else np.nan

    se_boot = se_ls
    if sigma is not None:
        s = np.asarray(sigma, dtype=np.float64)[m]
        s = np.where(np.isfinite(s) & (s > 0), s, np.nan)
        if np.isfinite(s).all():
            rng = np.random.default_rng(seed)
            sl = np.empty(n_mc, dtype=np.float64)
            for i in range(n_mc):
                yy = y + rng.standard_normal(y.size) * s
                sl[i] = ((x - xm) * (yy - yy.mean())).sum() / sxx
            se_boot = float(sl.std(ddof=1))

    t = b / se_ls if np.isfinite(se_ls) and se_ls > 0 else np.nan
    return dict(n=int(L.size), slope=b, se_ls=se_ls, se_boot=se_boot, t=t,
                significant=bool(np.isfinite(t) and b > 0 and t > SIGMA_SIG))


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="fgnspec")
    ap.add_argument("--rule", default="bd")
    ap.add_argument("--checkpoints", type=int, default=4)
    ap.add_argument("--if-boot", type=int, default=16,
                    help="random reference-run resplits of the Isolation Forest per cell "
                         "and variant (analyze.py default: 16)")
    ap.add_argument("--if-trees", type=int, default=400,
                    help="n_estimators (analyze.py default: 400)")
    ap.add_argument("--if-runs", type=int, default=96,
                    help="subsample every cell to this many runs before the "
                         "Isolation Forest step; 0 = use all runs (analyze.py "
                         "default: 0, README documents 96 for the campaigns)")
    ap.add_argument("--seed", type=int, default=7, help="analyze.py default: 7")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    variants = [v.strip().upper() for v in args.variants.split(",") if v.strip()]
    for v in variants:
        if v not in VARIANTS:
            raise SystemExit(f"unknown variant {v!r}; use {VARIANTS}")

    files = discover(args.tag, args.rule)
    archived = read_diagnostics(args.tag, args.rule)
    cells = sorted(files.keys(), key=lambda k: (k[1], k[0]))
    thetas = sorted({k[0] for k in cells})
    Ls = sorted({k[1] for k in cells})
    print(f"tag={args.tag} rule={args.rule}  {len(cells)} raw cells  "
          f"thetas={thetas}  Ls={Ls}")
    print(f"archived diagnostics: {len(archived)} rows")
    print(f"variants={variants}  checkpoints={args.checkpoints} "
          f"if_runs={args.if_runs} if_boot={args.if_boot} "
          f"if_trees={args.if_trees} seed={args.seed}")
    print(f"removed in NO_INC     : {sorted(INCREMENT_FEATURES)}")
    print("removed in NO_SCALING : alpha_loc_r*, logCnorm_r*, alpha_s")
    print("kept in every variant : skew_h, kurt_h")

    # ---- one feature pass per cell, then the three column subsets ----
    feats = {}          # (theta, L) -> (names, F)   F is (n_runs, n_ck, n_f)
    for (theta, L) in cells:
        t0 = time.time()
        res = load_result(files[(theta, L)])
        nc = res["profiles"].shape[1]
        k0 = max(0, nc - args.checkpoints)
        names, F, _ = snapshot_features(res["profiles"][:, k0:, :])
        feats[(theta, L)] = (names, F)
        print(f"  features theta={theta:+.4f} L={L:5d}  "
              f"profiles={res['profiles'].shape}  F={F.shape}  n_f={len(names)}"
              f"  ({time.time() - t0:.1f}s)", flush=True)

    # ---- Isolation Forest for every variant ----
    fit = {}            # (theta, L, variant) -> dict
    for L in Ls:
        if (0.0, L) not in feats:
            print(f"  !! no theta=0 reference for L={L}; L skipped (as in "
                  f"analyze.py)")
            continue
        names_ref, F_ref = feats[(0.0, L)]
        for theta in thetas:
            if (theta, L) not in feats:
                continue
            names_th, F_th = feats[(theta, L)]
            if names_th != names_ref:
                raise SystemExit(f"feature names differ for theta={theta} "
                                 f"L={L}")
            nr_th = F_th.shape[0]
            seed = args.seed + int(round(theta * 1e4)) + 10 ** 6

            m = min(args.if_runs or nr_th, nr_th, F_ref.shape[0])
            for var in variants:
                cols = variant_columns(names_ref, var)
                # Keep run identity; the detector averages late snapshots per run.
                Xr = F_ref[:m][:, :, cols]
                Xt = Xr if theta == 0 else F_th[:m][:, :, cols]
                t0 = time.time()
                r = iforest_excess(Xr, Xt, args.if_trees, seed=seed,
                                   n_boot=args.if_boot)
                fit[(theta, L, var)] = dict(
                    zPhi=r["z_anom"], zPhi_err=r["z_anom_err"],
                    dPhi=r["dPhi"], dPhi_err=r["dPhi_err"],
                    n_features=int(cols.size), n_rows=int(Xr.shape[0]),
                    n_train=r["n_train"], n_calibration=r["n_calibration"], n_null=r["n_null"],
                    calibration_resolution=r["calibration_resolution"])
                print(f"  th={theta:+.4f} L={L:5d} {var:>10s} "
                      f"n_f={cols.size:2d} rows={Xr.shape[0]:4d} "
                      f"zPhi={r['z_anom']:+8.3f}+-{r['z_anom_err']:.3f} "
                      f"dPhi={r['dPhi']:+.4f}  ({time.time() - t0:.1f}s)",
                      flush=True)

    # ---- per-cell table over the cells present in BOTH diagnostics ----
    keys = sorted({(t, L) for (t, L, v) in fit if (t, L) in archived},
                  key=lambda k: (k[1], k[0]))
    print(f"\ncommon cells (raw features + archived Hellinger): {len(keys)}")

    def col(var, attr="zPhi"):
        return np.array([fit[(t, L, var)][attr] if (t, L, var) in fit else np.nan
                         for (t, L) in keys], dtype=np.float64)

    th = np.array([t for (t, L) in keys], dtype=np.float64)
    LL = np.array([L for (t, L) in keys], dtype=np.float64)
    zH = np.array([archived[k]["zH"] for k in keys], dtype=np.float64)
    zPhi_arch = np.array([archived[k]["zPhi_archived"] for k in keys])
    zPhi_arch_err = np.array([archived[k]["zPhi_archived_err"] for k in keys])
    zP = {v: col(v) for v in variants}
    zPe = {v: col(v, "zPhi_err") for v in variants}

    outdir = args.outdir or os.path.join(ROOT, "data", "analysis")
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, f"ablation_{args.tag}.csv")
    json_path = os.path.join(outdir, f"ablation_{args.tag}.json")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["theta", "L", "zH", "zPhi_FULL", "zPhi_NO_INC",
                    "zPhi_NO_SCALING"])
        for i, (t, L) in enumerate(keys):
            w.writerow([f"{t:.4f}", L, repr(float(zH[i]))]
                       + [repr(float(zP[v][i])) if v in zP else "nan"
                          for v in VARIANTS])

    # ---- correlations ----
    def corr_block(sel):
        out = {}
        for v in variants:
            out[v] = corr_pair(zH[sel], zP[v][sel])
        out["ARCHIVED"] = corr_pair(zH[sel], zPhi_arch[sel])
        return out

    all_sel = np.ones(len(keys), dtype=bool)
    pos_sel = th > 0
    corr_all = corr_block(all_sel)
    corr_pos = corr_block(pos_sel)

    print("\n" + "=" * 88)
    print("CORRELATION BETWEEN zH (Hellinger) AND zPhi (Isolation Forest)")
    print("=" * 88)
    for label, blk, nsel in (("all common cells", corr_all, int(all_sel.sum())),
                             ("theta > 0 only", corr_pos, int(pos_sel.sum()))):
        print(f"\n  [{label}]  n={nsel}")
        print(f"    {'zPhi source':>10s} {'Pearson':>10s} {'p':>10s} "
              f"{'Spearman':>10s} {'p':>10s}")
        for v in list(variants) + ["ARCHIVED"]:
            c = blk[v]
            print(f"    {v:>10s} {c['pearson']:+10.4f} {c['pearson_p']:10.2e} "
                  f"{c['spearman']:+10.4f} {c['spearman_p']:10.2e}")

    # ---- FULL vs ablated changes ----
    print("\n" + "=" * 88)
    print("CHANGE IN zPhi INDUCED BY THE ABLATION (paired, same n_boot/seed)")
    print("=" * 88)
    changes = {}
    print(f"  {'comparison':>26s} {'mean dz':>10s} {'RMS dz':>10s} "
          f"{'max |dz|':>10s} {'mean zPhi err':>14s}")
    for v in variants:
        if v == "FULL":
            continue
        d = zP[v] - zP["FULL"]
        d = d[np.isfinite(d)]
        changes[v] = dict(mean=float(d.mean()),
                          rms=float(np.sqrt((d ** 2).mean())),
                          max_abs=float(np.abs(d).max()), n=int(d.size),
                          mean_zphi_err=float(np.nanmean(zPe[v])))
        print(f"  {'FULL -> ' + v:>26s} {d.mean():+10.4f} "
              f"{np.sqrt((d ** 2).mean()):10.4f} {np.abs(d).max():10.4f} "
              f"{np.nanmean(zPe[v]):14.4f}")
    sd_full = float(np.nanstd(zP["FULL"]))
    print(f"\n  context: sd(zPhi_FULL) over the {len(keys)} cells = {sd_full:.4f};"
          f"  mean zPhi_FULL = {np.nanmean(zP['FULL']):+.4f};"
          f"  mean |zPhi_FULL| = {np.nanmean(np.abs(zP['FULL'])):.4f}")

    # ---- fixed-theta scale dependence ----
    print("\n" + "=" * 88)
    print("FIXED-THETA SCALE DEPENDENCE:  slope of z vs ln L  (least squares)")
    print("=" * 88)
    slopes = {}
    series = {"zH": (zH, None)}
    for v in variants:
        series["zPhi_" + v] = (zP[v], zPe[v])
    series["zPhi_ARCHIVED"] = (zPhi_arch, zPhi_arch_err)

    for name, (y, sig) in series.items():
        tab = {}
        for t in thetas:
            s = th == t
            if s.sum() < 3:
                continue
            tab[f"{t:+.4f}"] = slope_vs_lnL(
                LL[s], y[s], None if sig is None else sig[s],
                seed=int(round(t * 1e4)) + 3)
        slopes[name] = tab

    order = ["zH"] + ["zPhi_" + v for v in variants] + ["zPhi_ARCHIVED"]
    print("\n  slope d z / d ln L   (+- = least-squares sigma; * = >2 sigma, >0)")
    print(f"  {'series':>16s}" + "".join(f"{t:>16.4f}" for t in REPORT_THETAS))
    for name in order:
        line = []
        for t in REPORT_THETAS:
            r = slopes[name].get(f"{t:+.4f}")
            if r is None or not np.isfinite(r["slope"]):
                line.append(f"{'--':>16s}")
            else:
                mark = "*" if r["significant"] else " "
                line.append(f"{r['slope']:+9.3f}+-{r['se_ls']:.3f}{mark}")
        print(f"  {name:>16s}" + "".join(line))

    print("\n  full theta grid, all series"
          "  (slope +- se_ls, se_boot, t, significant_positive)")
    for name in order:
        print(f"\n    {name}")
        for t in thetas:
            r = slopes[name].get(f"{t:+.4f}")
            if r is None:
                continue
            print(f"      theta={t:+.4f}  n={r['n']}  slope={r['slope']:+8.4f}"
                  f" +- {r['se_ls']:.4f} (ls)  +- {r['se_boot']:.4f} (boot)"
                  f"  t={r['t']:+6.2f}  significant+={r['significant']}")

    # ---- qualitative pattern comparison ----
    def pattern(name):
        return tuple(t for t in thetas
                     if slopes[name].get(f"{t:+.4f}", {}).get("significant"))

    print("\n" + "=" * 88)
    print("QUALITATIVE PATTERN: theta values with a significant (t > 2) positive"
          " slope of z against ln L")
    print("=" * 88)
    pats = {name: pattern(name) for name in slopes}
    for name in order:
        pl = ", ".join(f"{t:.4f}" for t in pats[name]) or "(none)"
        print(f"  {name:>16s}: {pl}")
    same_inc = pats["zPhi_FULL"] == pats["zPhi_NO_INC"]
    same_sc = pats["zPhi_FULL"] == pats["zPhi_NO_SCALING"]
    same_all = same_inc and same_sc
    print(f"\n  NO_INC     reproduces the FULL pattern exactly: {same_inc}")
    print(f"  NO_SCALING reproduces the FULL pattern exactly: {same_sc}")
    print(f"  => qualitative scale-dependence pattern unchanged by the ablation:"
          f" {same_all}")
    if not same_all:
        for name in order:
            only = set(pats[name]) - set(pats["zPhi_FULL"])
            missing = set(pats["zPhi_FULL"]) - set(pats[name])
            print(f"     {name:>16s}: extra={sorted(only) or 'none'}  "
                  f"lost={sorted(missing) or 'none'}")

    # ---- json ----
    summary = dict(
        tag=args.tag, rule=args.rule, script="code/feature_ablation.py",
        config=dict(checkpoints=args.checkpoints, if_boot=args.if_boot,
                    if_trees=args.if_trees, if_runs=args.if_runs,
                    seed=args.seed, alpha=0.02, variants=variants,
                    only_difference="which feature columns are retained"),
        statistics=dict(version="independent-runs-v1", unit="run",
                        features="late-snapshot mean within run",
                        null="reference runs held out of training and threshold calibration",
                        errors="dispersion across independent-run reference resplits",
                        standardized_scores_are_significance=False,
                        scale_patterns_are_physical_crossovers=False),
        feature_families=dict(
            increment=sorted(INCREMENT_FEATURES),
            scaling=["alpha_loc_r*", "logCnorm_r*", "alpha_s"],
            kept_in_every_variant=["skew_h", "kurt_h"]),
        n_features_by_L={str(L): {v: int(variant_columns(
            feats[(0.0, L)][0], v).size) for v in variants}
            for L in Ls if (0.0, L) in feats},
        cells=[dict(theta=float(t), L=int(L), zH=float(zH[i]),
                    zPhi_archived=float(zPhi_arch[i]),
                    zPhi_archived_err=float(zPhi_arch_err[i]),
                    **{f"zPhi_{v}": float(zP[v][i]) for v in variants},
                    **{f"zPhi_{v}_err": float(zPe[v][i]) for v in variants},
                    **{f"dPhi_{v}": float(fit[(t, L, v)]["dPhi"])
                       for v in variants},
                    **{k: fit[(t, L, variants[0])][k] for k in
                       ("n_rows", "n_train", "n_calibration", "n_null", "calibration_resolution")})
               for i, (t, L) in enumerate(keys)],
        correlations_all_cells=corr_all,
        correlations_theta_positive=corr_pos,
        ablation_changes=changes,
        sd_zPhi_FULL=sd_full,
        slopes_vs_lnL=slopes,
        significant_positive_slope_thetas={k: list(v) for k, v in pats.items()},
        verdict=dict(pattern_unchanged_by_NO_INC=same_inc,
                     pattern_unchanged_by_NO_SCALING=same_sc,
                     pattern_unchanged=same_all),
    )
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nwrote {csv_path}")
    print(f"wrote {json_path}")

    # ---- compact final summary ----
    print("\n" + "#" * 88)
    print("# COMPACT SUMMARY")
    print("#" * 88)
    print(f"  cells: {len(keys)}  (theta {thetas[0]:.4f}..{thetas[-1]:.4f}, "
          f"L {Ls[0]}..{Ls[-1]})   n_boot={args.if_boot} "
          f"n_estimators={args.if_trees} if_runs={args.if_runs} "
          f"seed={args.seed} checkpoints={args.checkpoints}")
    for label, blk in (("all", corr_all), ("theta>0", corr_pos)):
        print("  corr(zH, zPhi) [%s, n=%d]: " % (label, blk[variants[0]]["n"])
              + "  ".join(f"{v}: r={blk[v]['pearson']:+.3f} "
                          f"rho={blk[v]['spearman']:+.3f}" for v in variants))
    for v in variants:
        if v == "FULL":
            continue
        c = changes[v]
        print(f"  FULL->{v:>10s}: mean dz={c['mean']:+.4f}  "
              f"RMS dz={c['rms']:.4f}  (per-cell zPhi err ~"
              f"{c['mean_zphi_err']:.3f})")
    print(f"  significant positive slopes: FULL={pats['zPhi_FULL']} "
          f"NO_INC={pats['zPhi_NO_INC']} NO_SCALING={pats['zPhi_NO_SCALING']}"
          f"  zH={pats['zH']}")
    print(f"  pattern unchanged by ablation: {same_all}")


if __name__ == "__main__":
    main()
