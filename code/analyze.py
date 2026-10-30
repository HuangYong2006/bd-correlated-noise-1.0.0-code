"""
Analysis driver.

Produces, for each system size L:
  * dH(theta,L)  -- bias-corrected Hellinger distance between exact
                    height-INCREMENT pmfs of the theta ensemble and the
                    theta = 0 (KPZ) reference ensemble   [diagnostic 1, primary]
  * dH_h         -- the same for the height density itself [control; expected
                    to be blind, which is itself a reported result]
  * zPhi         -- standardised mean Isolation-Forest anomaly score of the
                    amplitude-free scaling/morphology features w.r.t. the
                    theta = 0 reference                        [diagnostic 2]
  * Phi          -- the calibrated anomalous FRACTION (secondary iForest output)
  * alpha_s, alpha_loc(r), slope/facet descriptors, W_sat, beta
Legacy threshold crossings and fits are retained as descriptive summaries.
Their pointwise uncertainty approximation ignores shared-reference covariance
and must not be interpreted as physical transition estimates.

Usage:  python code/analyze.py --tag main --rule bd --checkpoints 4
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from bd_core import load_result
from features import snapshot_features, alpha_loc_at
from diagnostics import (hellinger_excess, iforest_excess, onset_theta,
                         fit_onset_scaling, hellinger_null_curve,
                         increment_pmf, terrace_length_pmf,
                         discrete_hellinger_excess, onset_with_error)


# --------------------------------------------------------------------------

def discover(tag, rule):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if not m:
            continue
        tok = m.group("tok")
        sign = -1.0 if tok[0] == "m" else 1.0
        out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = f
    return out


def beta_in_window(res, frac_lo, frac_hi):
    """Growth exponent in t in [frac_lo, frac_hi] * L^1.5 (t_x ~ 0.6 L^1.5)."""
    L = res["spec"]["L"]
    wt = res["w_times"].astype(np.float64)
    W = np.sqrt(res["w_series"].astype(np.float64).mean(axis=0))
    lo, hi = frac_lo * L ** 1.5, frac_hi * L ** 1.5
    m = (wt >= lo) & (wt <= hi) & (W > 0)
    if m.sum() < 5:
        return np.nan
    return float(np.polyfit(np.log(wt[m]), np.log(W[m]), 1)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--rule", default="bd")
    ap.add_argument("--checkpoints", type=int, default=4)
    ap.add_argument("--bins", type=int, default=64)
    ap.add_argument("--splits", type=int, default=120)
    ap.add_argument("--hell-boot", type=int, default=200,
                    help="outer independent-run bootstraps of the complete exact-pmf excess")
    ap.add_argument("--if-boot", type=int, default=16,
                    help="random disjoint reference-run resplits (not bootstrap standard errors)")
    ap.add_argument("--if-trees", type=int, default=400)
    ap.add_argument("--if-runs", type=int, default=0,
                    help="subsample every cell to this many runs before the "
                         "Isolation Forest step, so that the number of TRAINING "
                         "samples is identical across system sizes and across "
                         "theta.  Without this the detector's statistical power "
                         "varies with L, and any significance-based onset "
                         "criterion drifts with ensemble size rather than with "
                         "physics.  0 = use all runs (not recommended for "
                         "finite-size scaling).")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--r-frac", type=float, default=0.125,
                    help="probe separation for the Hellinger 'diff' observable, "
                         "as a fraction of L (keeps the probed scale comparable "
                         "across system sizes for the finite-size-scaling test)")
    ap.add_argument("--hell-runs", type=int, default=0,
                    help="subsample every cell to this many runs before the "
                         "Hellinger step.  The split-half null floor of the "
                         "distributional diagnostic scales with the number of "
                         "runs, so a mismatch across L biases the extracted "
                         "onset with statistical power rather than physics.  "
                         "0 = use all runs.")
    ap.add_argument("--skip-iforest", action="store_true",
                    help="skip the (expensive) Isolation Forest step; useful for "
                         "fast iteration on the distributional diagnostic, the "
                         "onsets and the exponents")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    files = discover(args.tag, args.rule)
    if not files:
        print("no cells found"); return
    cells = {k: load_result(v) for k, v in files.items()}
    thetas = sorted({k[0] for k in cells})
    Ls = sorted({k[1] for k in cells})
    print(f"tag={args.tag} rule={args.rule}  {len(cells)} cells  thetas={thetas}  Ls={Ls}")

    outdir = args.outdir or os.path.join(ROOT, "data", "analysis", f"{args.tag}_{args.rule}")
    os.makedirs(outdir, exist_ok=True)

    rows = []
    for L in Ls:
        if (0.0, L) not in cells:
            print(f"  !! no theta=0 reference for L={L}"); continue
        ref = cells[(0.0, L)]
        nc = ref["profiles"].shape[1]
        k0 = max(0, nc - args.checkpoints)
        def _sub(arr):
            if not args.hell_runs:
                return arr
            return arr[:min(args.hell_runs, arr.shape[0])]

        ref_last = _sub(ref["profiles"][:, -1, :])
        names, F_ref, extra_ref = snapshot_features(ref["profiles"][:, k0:, :])
        r_probe = max(2, int(args.r_frac * L))
        # --- primary distributional observable: exact, binning-free pmfs ---
        ref_inc = increment_pmf(ref_last)
        ref_terr = terrace_length_pmf(ref_last)
        n_i = discrete_hellinger_excess(ref_inc, ref_inc, args.splits, args.seed, n_boot=0)
        n_t = discrete_hellinger_excess(ref_terr, ref_terr, args.splits, args.seed + 1, n_boot=0)
        # --- controls: histogrammed variants (documented as unsuitable) ---
        hn_d, hn_dsd = hellinger_null_curve(ref_last, args.bins, args.splits,
                                            args.seed, "diff", r_probe)
        hn_h, hn_hsd = hellinger_null_curve(ref_last, args.bins, args.splits,
                                            args.seed, "height", 1)
        print(f"  [L={L}] null floor: incr_pmf={n_i['Hnull_mean']:.5f}"
              f"+-{n_i['Hnull_std']:.5f}  terr_pmf={n_t['Hnull_mean']:.5f}"
              f"+-{n_t['Hnull_std']:.5f}  | hist diff(r={r_probe})={hn_d:.4f}"
              f"  hist height={hn_h:.4f}", flush=True)

        for theta in thetas:
            if (theta, L) not in cells:
                continue
            r = cells[(theta, L)]
            th_last = _sub(r["profiles"][:, -1, :])
            hd = discrete_hellinger_excess(ref_inc, ref_inc if theta == 0 else increment_pmf(th_last),
                                           args.splits, args.seed, n_boot=args.hell_boot)
            ht = discrete_hellinger_excess(ref_terr, ref_terr if theta == 0 else terrace_length_pmf(th_last),
                                           args.splits, args.seed + 1, n_boot=args.hell_boot)
            # Histogram controls get conditional null floors only; no inferential
            # error is claimed for these documented unsuitable observables.
            hx = hellinger_excess(ref_last, th_last, args.bins, args.splits,
                                  args.seed, None, "diff", r_probe, n_boot=0)
            hh = hellinger_excess(ref_last, th_last, args.bins, args.splits,
                                  args.seed, None, "height", 1, n_boot=0)
            _, F_th, extra_th = snapshot_features(r["profiles"][:, k0:, :])
            nr_th = F_th.shape[0]
            m = min(args.if_runs or nr_th, nr_th, F_ref.shape[0])
            # Subsample runs first. iforest_excess averages snapshots within run.
            Xr = F_ref[:m]
            X_th = Xr if theta == 0 else F_th[:m]
            ires = iforest_excess(Xr, X_th, args.if_trees,
                                  seed=args.seed + int(round(theta * 1e4)) + 10 ** 6,
                                  n_boot=args.if_boot) if not args.skip_iforest else \
                dict(z_anom=np.nan, z_anom_err=np.nan, dPhi=np.nan, dPhi_err=np.nan,
                     n_ref_runs=m, n_theta_runs=m, n_train=0, n_calibration=0, n_null=0,
                     calibration_resolution=np.nan, phi_null=np.nan,
                     phi_null_ci_low=np.nan, phi_null_ci_high=np.nan)
            W_sat = float(np.sqrt(r["w_series"].astype(np.float64).mean(axis=0))[-1])

            rows.append(dict(
                theta=theta, L=L,
                dH=hd["dH"], dH_err=hd["dH_err"], zH=hd["z"], zH_err=np.nan,
                dH_terr=ht["dH"], dH_terr_err=ht["dH_err"], zH_terr=ht["z"],
                zH_terr_err=np.nan,
                dH_histdiff=hx["dH"], zH_histdiff=hx["z"], zH_histdiff_err=np.nan,
                dH_height=hh["dH"], zH_height=hh["z"], zH_height_err=np.nan,
                Hsig_split_std=hd["Hsig_split_std"], dH_split_std=hd["dH_split_std"],
                hell_n_ref_runs=hd["n_ref_runs"], hell_n_theta_runs=hd["n_theta_runs"],
                hell_n_bootstrap=hd["n_bootstrap"],
                Hnull=n_i["Hnull_mean"], Hnull_std=n_i["Hnull_std"],
                Hnull_terr=n_t["Hnull_mean"], Hnull_histdiff=hn_d,
                Hnull_height=hn_h, r_probe=r_probe,
                zPhi=ires["z_anom"], zPhi_err=ires["z_anom_err"],
                dPhi=ires["dPhi"], dPhi_err=ires["dPhi_err"],
                if_n_ref_runs=ires["n_ref_runs"], if_n_theta_runs=ires["n_theta_runs"],
                if_n_train=ires["n_train"], if_n_calibration=ires["n_calibration"],
                if_n_null=ires["n_null"], if_calibration_resolution=ires["calibration_resolution"],
                if_phi_null=ires["phi_null"], if_phi_null_ci_low=ires["phi_null_ci_low"],
                if_phi_null_ci_high=ires["phi_null_ci_high"],
                W_sat=W_sat,
                alpha_s=float(np.nanmean(extra_th["alpha_s"])),
                alpha_loc_r8=alpha_loc_at(extra_th, 8),
                alpha_loc_r32=alpha_loc_at(extra_th, 32),
                alpha_loc_r128=alpha_loc_at(extra_th, 128),
                skew_d=float(np.mean(extra_th["skew_d"])),
                kurt_d=float(np.mean(extra_th["kurt_d"])),
                frac_flat=float(np.mean(extra_th["frac_flat"])),
                frac_step2=float(np.mean(extra_th["frac_step2"])),
                frac_step3=float(np.mean(extra_th["frac_step3"])),
                mean_flat_run=float(np.mean(extra_th["mean_flat_run"])),
                log_mean_flat_run=float(np.mean(np.log10(np.maximum(extra_th["mean_flat_run"], 1e-6)))),
                slope_entropy=float(np.mean(extra_th["slope_entropy"])),
                skew_h=float(np.mean(extra_th["skew_h"])),
                kurt_h=float(np.mean(extra_th["kurt_h"])),
                beta_1=beta_in_window(r, 0.01, 0.10),
                beta_2=beta_in_window(r, 0.03, 0.25),
                beta_3=beta_in_window(r, 0.10, 0.60),
            ))
            print(f"  th={theta:+.4f} L={L:5d} dHpmf={hd['dH']:+.5f} zH={hd['z']:7.2f} "
                  f"| terr={ht['dH']:+.5f} zT={ht['z']:6.2f} "
                  f"| histdiff={hx['dH']:+.4f} height={hh['dH']:+.4f} "
                  f"| zPhi={ires['z_anom']:8.2f} Phi={ires['dPhi']:+.3f} "
                  f"| W={W_sat:7.2f} as={np.mean(extra_th['alpha_s']):.3f} "
                  f"kd={np.mean(extra_th['kurt_d']):.2f}",
                  flush=True)

    # ---------------- onset + finite-size scaling ----------------
    print("\n" + "=" * 106)
    print("DESCRIPTIVE CROSSINGS AND FITS (not physical transition estimates)")
    print("=" * 106)
    summary = {"tag": args.tag, "rule": args.rule, "Ls": Ls, "features": names,
               "statistics": dict(version="independent-runs-v1", unit="run",
                   forest_features="late-snapshot mean within run",
                   forest_error="reference-run resplit dispersion, not a sampling standard error",
                   forest_null="independent held-out reference runs",
                   forest_null_interval="envelope of conditional 95% exact binomial intervals",
                   hellinger_error="outer independent-run bootstrap of full excess",
                   histogram_controls="conditional null only; error not estimated",
                   onset_error="pointwise Gaussian approximation; shared-reference covariance omitted",
                   onsets_are_physical=False)}
    onset_tab = {}
    criteria = [
        # (label, column, thresholds)
        # effect-size criteria: statistics-independent, hence suitable for a
        # drift analysis; significance (z) criteria are reported as well but
        # their location moves with the ensemble size by construction.
        ("Hellinger_pmf_effectsize", "dH", [0.005, 0.01, 0.02]),
        ("Hellinger_pmf_z", "zH", [2.0, 3.0, 5.0, 10.0]),
        ("Hellinger_terrace_z", "zH_terr", [2.0, 3.0, 5.0]),
        ("iForest_frac_effectsize", "dPhi", [0.05, 0.10, 0.20, 0.50]),
        ("iForest_zscore", "zPhi", [3.0, 5.0, 10.0]),
        ("Hellinger_histdiff_z", "zH_histdiff", [3.0, 5.0]),
        ("Hellinger_height_z", "zH_height", [3.0, 5.0, 10.0]),
        ("frac_step2", "frac_step2", [0.01, 0.05, 0.10]),
    ]
    for obs, key, thrs in criteria:
        onset_tab[obs] = {}
        for zthr in thrs:
            th_star, th_err = [], []
            for L in Ls:
                sub = sorted([r for r in rows if r["L"] == L], key=lambda r: r["theta"])
                th = np.array([r["theta"] for r in sub])
                val = np.array([r[key] for r in sub], dtype=np.float64)
                # per-point uncertainty: use the calibrated null scatter where
                # available, else the bootstrap spread of the observable
                sig = np.array([r.get(key + "_err", np.nan) for r in sub],
                               dtype=np.float64)
                ob = onset_with_error(th, val, sig, zthr, True, 300, 11)
                th_star.append(ob["theta_star"])
                th_err.append(ob["theta_star_err"])
            th_err_a = np.array([np.nan if e is None else e for e in th_err],
                                dtype=np.float64)
            fit = fit_onset_scaling(Ls, th_star, th_err_a, theta_c_free=True)
            fit0 = fit_onset_scaling(Ls, th_star, th_err_a, theta_c_free=False)
            onset_tab[obs][f"thr{zthr}"] = dict(                theta_star=[None if not np.isfinite(t) else float(t) for t in th_star],
                free={k: (None if v is None or (isinstance(v, float) and not np.isfinite(v)) else v)
                      for k, v in fit.items()},
                forced_zero={k: (None if v is None or (isinstance(v, float) and not np.isfinite(v)) else v)
                             for k, v in fit0.items()})
            s1 = ", ".join(f"{L}:{t:.4f}" if np.isfinite(t) else f"{L}:--"
                           for L, t in zip(Ls, th_star))
            e1 = ", ".join(f"+-{e:.4f}" if e is not None and np.isfinite(e) else ""
                           for e in th_err)
            print(f"\n  [{obs}] criterion {zthr}:  theta*(L) = {s1}")
            print(f"      onset errors: {e1}")
            print(f"      free   theta_c = {fit['theta_c']:.4f} +- {fit['theta_c_err']:.4f}"
                  f"   1/nu = {1/fit['nu']:.3f}   chi2/dof = {fit['chi2']:.2f}/{fit['dof']}"
                  f"   p = {fit['pvalue']:.3f}")
            print(f"      th_c=0 forced     A = {fit0['A']:.4f}"
                  f"   1/nu = {1/fit0['nu']:.3f}   chi2/dof = {fit0['chi2']:.2f}/{fit0['dof']}"
                  f"   p = {fit0['pvalue']:.3f}")
    summary["onsets"] = onset_tab

    # ---------------- direct exponents ----------------
    print("\n" + "=" * 106)
    print("DIRECT EXPONENTS  (alpha from W_sat ~ L^alpha, beta in growth window)")
    print("=" * 106)
    exp_tab = {}
    for theta in thetas:
        Ws, Ls_ok = [], []
        for L in Ls:
            sub = [r for r in rows if r["L"] == L and r["theta"] == theta]
            if sub:
                Ws.append(sub[0]["W_sat"]); Ls_ok.append(L)
        if len(Ls_ok) >= 3:
            p = np.polyfit(np.log(Ls_ok), np.log(Ws), 1)
            resid = np.log(Ws) - np.polyval(p, np.log(Ls_ok))
            alpha = float(p[0])
            alpha_err = float(np.sqrt((resid ** 2).sum() / max(1, len(Ls_ok) - 2)
                                      / max(1e-9, ((np.log(Ls_ok) - np.log(Ls_ok).mean()) ** 2).sum())))
        else:
            alpha, alpha_err = np.nan, np.nan
        b2 = [r["beta_2"] for r in rows if r["theta"] == theta]
        b1 = [r["beta_1"] for r in rows if r["theta"] == theta]
        exp_tab[f"{theta:+.4f}"] = dict(
            alpha=alpha, alpha_err=alpha_err,
            alpha_error_method="log-fit residual error, not run-bootstrap uncertainty",
            beta=float(np.nanmean(b2)), beta_early=float(np.nanmean(b1)),
            W_by_L={str(L): float(r["W_sat"]) for r in rows if r["theta"] == theta})
        print(f"  theta={theta:+.4f}  alpha={alpha:6.4f}+-{alpha_err:.4f}"
              f"   beta={np.nanmean(b2):6.4f}   beta(early)={np.nanmean(b1):6.4f}")
    summary["exponents"] = exp_tab

    keys = list(rows[0].keys())
    with open(os.path.join(outdir, "diagnostics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    np.savez(os.path.join(outdir, "rows.npz"),
             **{k: np.array([r[k] for r in rows], dtype=np.float64) for k in keys})
    print(f"\nwrote {outdir}/diagnostics.csv, summary.json, rows.npz")


if __name__ == "__main__":
    main()
