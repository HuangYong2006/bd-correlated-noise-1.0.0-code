"""
Figure generation for the manuscript.

Reads data/raw/<tag>/ (simulation cells) and data/analysis/<tag>_<rule>/
(diagnostics tables) and writes vector PDF + PNG figures to figs/<tag>/.

Usage:  python code/figures.py --tag main --rule bd
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from bd_core import load_result, ffgn_weights, ffgn_n_components, empirical_correlation

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.linewidth": 1.3,
})

THETA_C_RG = 1.0 / 6.0
THETA_C_LIT = 0.25


def discover(tag, rule):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = f
    return out


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  figs/{os.path.basename(outdir)}/{name}.pdf/.png")


# --------------------------------------------------------------------------

def fig_noise(outdir):
    """Noise-generator characterisation."""
    from numba import njit
    from bd_core import _u01

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

    T, n_sites = 1 << 20, 12
    max_lag = 1 << 13
    thetas = [0.05, 0.10, 0.1667, 0.25, 0.30, 0.40]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))
    tau = np.arange(1, max_lag)
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(thetas)))
    rec = []
    for c, th in zip(cols, thetas):
        N = ffgn_n_components(T)
        Wn, rn = ffgn_weights(th, N)
        for _ in range(n_sites):
            pass
        v = ffgn_seq(T, Wn, rn, seed=11)
        # fast autocorrelation via FFT
        x = v - v.mean()
        nfft = 1 << int(np.ceil(np.log2(2 * T)))
        F = np.fft.rfft(x, n=nfft)
        ac = np.fft.irfft(F * np.conj(F), n=nfft)[:T]
        ac /= ac[0]
        axes[0].loglog(tau, np.abs(ac[1:max_lag]), color=c, lw=0.9,
                       label=rf"$\theta={th:.4f}$")
        m = tau <= 4096
        sl = np.polyfit(np.log(tau[m]), np.log(np.abs(ac[1:max_lag][m])), 1)[0]
        rec.append((th, 2 * th - 1, sl, (1 + sl) / 2))
    axes[0].set_xlabel(r"lag $\tau$  [sweeps]")
    axes[0].set_ylabel(r"$C_\eta(\tau)/C_\eta(0)$")
    axes[0].set_title("(a) FFGN temporal correlation")
    axes[0].legend(ncol=2, loc="lower left")

    rec = np.array(rec)
    axes[1].plot([0, 0.45], [0, 0.45], "k--", lw=1, label="nominal = effective")
    axes[1].errorbar(rec[:, 0], rec[:, 3], fmt="o", ms=4, color="C3",
                     label="measured $\\theta_{\\rm eff}$")
    axes[1].set_xlabel(r"nominal $\theta$")
    axes[1].set_ylabel(r"effective $\theta$ from $C_\eta$ slope")
    axes[1].set_title("(b) generator calibration")
    axes[1].legend(loc="upper left")

    # fraction of 1s in the binarised noise vs theta
    ths = np.linspace(0.0, 0.45, 10)
    fr = []
    for th in ths:
        if th == 0:
            rng = np.random.default_rng(3)
            fr.append((rng.random(1 << 16) > 0.5).mean())
        else:
            N = ffgn_n_components(1 << 16)
            Wn, rn = ffgn_weights(th, N)
            v = ffgn_seq(1 << 16, Wn, rn, seed=5)
            fr.append(float((v > 0).mean()))
    axes[2].plot(ths, fr, "o-", color="C0")
    axes[2].axhline(0.5, color="k", ls="--", lw=1)
    axes[2].set_ylim(0.45, 0.55)
    axes[2].set_xlabel(r"$\theta$")
    axes[2].set_ylabel(r"$P(\eta'=1)$")
    axes[2].set_title("(c) binarised noise has zero mean")
    fig.tight_layout()
    save(fig, outdir, "fig1_noise")
    return rec


def fig_theta0(cells, outdir):
    """theta = 0 validation: growth and saturation."""
    Ls = sorted({L for (_, L) in cells})
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))
    Wsat, Ls_ok = [], []
    for L in Ls:
        r = cells.get((0.0, L))
        if r is None:
            continue
        wt = r["w_times"].astype(np.float64)
        W = np.sqrt(r["w_series"].astype(np.float64).mean(axis=0))
        axes[0].loglog(wt, W, label=f"$L={L}$", lw=1.0)
        # guide slope 1/3
        axes[0].loglog(wt, 0.30 * wt ** (1 / 3), ls=":", color="k", lw=0.8)
        Wsat.append(W[-1])
        Ls_ok.append(L)
    axes[0].set_xlabel(r"$t$  [sweeps]")
    axes[0].set_ylabel(r"$W(L,t)$")
    axes[0].set_title(r"(a) growth: dotted line $\propto t^{1/3}$")
    axes[0].legend(ncol=2, fontsize=7)

    Ls_ok = np.array(Ls_ok, dtype=float)
    Wsat = np.array(Wsat)
    axes[1].loglog(Ls_ok, Wsat, "o", color="C0", label="saturated $W$")
    if Ls_ok.size >= 3:
        p = np.polyfit(np.log(Ls_ok), np.log(Wsat), 1)
        axes[1].loglog(Ls_ok, np.exp(np.polyval(p, np.log(Ls_ok))), "-",
                       color="C3", label=rf"fit $\alpha={p[0]:.3f}$")
    axes[1].set_xlabel(r"$L$")
    axes[1].set_ylabel(r"$W_{\rm sat}$")
    axes[1].set_title(r"(b) saturation: expect $\alpha=1/2$")
    axes[1].legend()

    lo, hi = [], []
    for L in Ls:
        r = cells.get((0.0, L))
        if r is None:
            continue
        wt = r["w_times"].astype(np.float64)
        W = np.sqrt(r["w_series"].astype(np.float64).mean(axis=0))
        m = (wt >= 0.01 * L ** 1.5) & (wt <= 0.25 * L ** 1.5)
        if m.sum() > 5:
            lo.append(L)
            hi.append(np.polyfit(np.log(wt[m]), np.log(W[m]), 1)[0])
    axes[2].semilogx(lo, hi, "s-", color="C2")
    axes[2].axhline(1 / 3, color="k", ls="--", lw=1, label=r"$\beta=1/3$")
    axes[2].set_xlabel(r"$L$")
    axes[2].set_ylabel(r"$\beta_{\rm eff}$ in growth window")
    axes[2].set_title("(c) measured growth exponent")
    axes[2].legend()
    fig.tight_layout()
    save(fig, outdir, "fig2_theta0")


def fig_widths(cells, outdir, L=1024):
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2))
    ths = sorted({th for (th, LL) in cells if LL == L})
    cols = plt.cm.plasma(np.linspace(0, 0.9, max(len(ths), 2)))
    for c, th in zip(cols, ths):
        r = cells[(th, L)]
        wt = r["w_times"].astype(np.float64)
        W = np.sqrt(r["w_series"].astype(np.float64).mean(axis=0))
        axes[0].loglog(wt, W, color=c, lw=1.0)
        axes[1].loglog(wt, W / W[-1], color=c, lw=1.0,
                       label=rf"$\theta={th:.3f}$")
    axes[0].set_xlabel(r"$t$  [sweeps]"); axes[0].set_ylabel(r"$W(L,t)$")
    axes[0].set_title(rf"(a) $L={L}$: width growth")
    axes[1].set_xlabel(r"$t$  [sweeps]"); axes[1].set_ylabel(r"$W(t)/W_{\rm sat}$")
    axes[1].set_title("(b) normalised")
    axes[1].legend(ncol=2, fontsize=6)
    fig.tight_layout()
    save(fig, outdir, f"fig3_widths_L{L}")


def fig_diagnostics(rows, outdir):
    """The two diagnostics versus theta, one panel per L."""
    Ls = sorted(set(rows["L"].astype(int)))
    fig, axes = plt.subplots(2, len(Ls), figsize=(2.6 * len(Ls), 4.8),
                             sharex="col", squeeze=False)
    for j, L in enumerate(Ls):
        m = rows["L"].astype(int) == L
        th = rows["theta"][m]
        o = np.argsort(th)
        th = th[o]
        ax = axes[0, j]
        ax.plot(th, rows["dH"][m][o], "o-", ms=3.5, color="C0")
        ax.fill_between(th, rows["Hnull"][m][o] - 2 * rows["Hnull_std"][m][o],
                        rows["Hnull"][m][o] + 2 * rows["Hnull_std"][m][o],
                        color="grey", alpha=0.4, lw=0)
        ax.axvline(THETA_C_RG, color="C3", ls="--", lw=1)
        ax.axvline(THETA_C_LIT, color="C2", ls=":", lw=1)
        ax.set_title(rf"$L={L}$")
        if j == 0:
            ax.set_ylabel(r"Hellinger excess $\Delta H$ (increment)")
        ax = axes[1, j]
        ax.plot(th, rows["dPhi"][m][o], "s-", ms=3.5, color="C1")
        ax.axvline(THETA_C_RG, color="C3", ls="--", lw=1)
        ax.axvline(THETA_C_LIT, color="C2", ls=":", lw=1)
        ax.set_xlabel(r"$\theta$")
        if j == 0:
            ax.set_ylabel(r"iForest anomaly fraction $\Phi$")
    axes[0, 0].plot([], [], color="C3", ls="--", label=r"$\theta_{\rm th}=1/6$")
    axes[0, 0].plot([], [], color="C2", ls=":", label=r"$\theta_c=0.25$")
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout()
    save(fig, outdir, "fig4_diagnostics")


def fig_onset_drift(summary, outdir):
    """The decisive figure: does theta*(L) drift to zero or to a finite value?"""
    Ls = np.array(summary["Ls"], dtype=float)
    onset = summary["onsets"]
    panels = [
        ("Hellinger_slope_effectsize", "dH", ["thr0.1", "thr0.2"]),
        ("iForest_frac_effectsize", "dPhi", ["thr0.2", "thr0.5"]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    for ax, (obs, _, keys) in zip(axes, panels):
        for k, mk in zip(keys, ("o", "s")):
            d = onset.get(obs, {}).get(k)
            if not d:
                continue
            ts = np.array([np.nan if v is None else v for v in d["theta_star"]],
                          dtype=float)
            m = np.isfinite(ts)
            if m.sum() < 2:
                continue
            ax.semilogx(Ls[m], ts[m], mk + "-", ms=5, label=k)
            fz = d.get("forced_zero") or {}
            if fz.get("A") is not None:
                Lg = np.geomspace(Ls[m].min(), Ls[m].max() * 4, 60)
                ax.semilogx(Lg, fz["A"] * Lg ** (-1.0 / fz["nu"]), "--", lw=1,
                            color="grey")
            ff = d.get("free") or {}
            if ff.get("theta_c") is not None:
                ax.axhline(ff["theta_c"], color="C3", lw=1, ls=":")
        ax.axhline(THETA_C_RG, color="C3", ls="--", lw=1.2)
        ax.axhline(THETA_C_LIT, color="C2", ls="-.", lw=1.2)
        ax.axhline(0.0, color="k", lw=1)
        ax.set_xlabel(r"$L$")
        ax.set_ylabel(r"onset $\theta^*(L)$")
        ax.legend(fontsize=7)
    axes[0].set_title("(a) Hellinger diagnostic")
    axes[1].set_title("(b) Isolation Forest diagnostic")
    axes[0].plot([], [], color="C3", ls="--", label=""); axes[0].legend(fontsize=7)
    fig.tight_layout()
    save(fig, outdir, "fig5_onset_drift")


def fig_exponents(rows, outdir):
    """alpha(theta), beta(theta) and the analytic branches."""
    Ls = sorted(set(rows["L"].astype(int)))
    ths = np.array(sorted(set(np.round(rows["theta"], 4))))
    alpha, alpha_err, beta = [], [], []
    for th in ths:
        W, LL = [], []
        for L in Ls:
            m = (rows["L"].astype(int) == L) & (np.round(rows["theta"], 4) == th)
            if m.sum():
                W.append(rows["W_sat"][m][0]); LL.append(L)
        if len(LL) >= 3:
            la = np.log(LL); lw = np.log(W)
            p = np.polyfit(la, lw, 1)
            alpha.append(p[0]); alpha_err.append(np.nan)
        else:
            alpha.append(np.nan); alpha_err.append(np.nan)
        m = np.round(rows["theta"], 4) == th
        beta.append(np.nanmean(rows["beta_2"][m]))
    alpha = np.array(alpha); beta = np.array(beta)

    t = np.linspace(0.0, 0.45, 200)
    med_x = (1 + 4 * t) / (3 + 2 * t)
    med_z = 5.0 / (3 + 2 * t)
    med_b = med_x / med_z
    mm_x = (1 + 2 * t) / (2 + t)
    mm_z = 3.0 / (2 + t)
    mm_b = mm_x / mm_z
    sc_x = 1.43 * t + 0.26

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.3))
    axes[0].errorbar(ths, alpha, yerr=alpha_err, fmt="o", ms=4, color="C0",
                     label=r"measured $\alpha$")
    axes[0].plot(t, med_x, "-", color="C3", lw=1.1, label="Medina one-loop")
    axes[0].plot(t, mm_x, "--", color="C2", lw=1.1, label="Ma & Ma (Flory)")
    axes[0].plot(t, sc_x, ":", color="C4", lw=1.1, label="Squizzato-Canet NPRG")
    axes[0].axhline(0.5, color="k", lw=0.8)
    axes[0].axvline(THETA_C_RG, color="grey", ls="--", lw=0.9)
    axes[0].set_xlabel(r"$\theta$"); axes[0].set_ylabel(r"$\alpha$")
    axes[0].set_title(r"(a) roughness exponent from $W_{\rm sat}\sim L^\alpha$")
    axes[0].legend(fontsize=7)

    axes[1].plot(ths, beta, "s-", ms=4, color="C0", label=r"measured $\beta$")
    axes[1].plot(t, med_b, "-", color="C3", lw=1.1, label="Medina one-loop")
    axes[1].plot(t, mm_b, "--", color="C2", lw=1.1, label="Ma & Ma (Flory)")
    axes[1].axhline(1 / 3, color="k", lw=0.8)
    axes[1].axvline(THETA_C_RG, color="grey", ls="--", lw=0.9)
    axes[1].set_xlabel(r"$\theta$"); axes[1].set_ylabel(r"$\beta$")
    axes[1].set_title("(b) growth exponent")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    save(fig, outdir, "fig6_exponents")


def fig_morphology(rows, cells, outdir, L=1024):
    ths = np.array(sorted(set(np.round(rows["theta"], 4))))
    m0 = rows["L"].astype(int) == L
    def series(col):
        v = []
        for th in ths:
            m = m0 & (np.round(rows["theta"], 4) == th)
            v.append(rows[col][m][0] if m.sum() else np.nan)
        return np.array(v)

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.4))
    ax = axes[0, 0]
    ax.plot(ths, series("alpha_s"), "o-", ms=4, color="C0")
    ax.set_ylabel(r"$\alpha_s$"); ax.set_title(r"(a) spectral exponent")
    for a in (0.5, 1.0, 1.5):
        ax.axhline(a, color="grey", lw=0.6, ls=":")
    ax = axes[0, 1]
    ax.semilogy(ths, np.maximum(series("frac_step2"), 1e-6), "s-", ms=4, color="C1")
    ax.set_ylabel(r"fraction of $|\Delta h|\geq 2$"); ax.set_title("(b) overhang/facet fraction")
    ax = axes[0, 2]
    ax.plot(ths, series("mean_flat_run"), "^-", ms=4, color="C2")
    ax.set_ylabel("mean terrace length"); ax.set_title("(c) terrace length")
    ax = axes[1, 0]
    ax.plot(ths, series("skew_d"), "o-", ms=4, color="C3")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"skew$(\Delta h)$")
    ax.set_title("(d) increment skewness")
    ax = axes[1, 1]
    ax.plot(ths, series("kurt_d"), "s-", ms=4, color="C4")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"kurt$(\Delta h)$")
    ax.set_title("(e) increment kurtosis")
    ax = axes[1, 2]
    for th in [0.0, 0.10, THETA_C_RG, 0.25, 0.40]:
        key = (round(th, 4), L)
        if key not in cells:
            continue
        h = cells[key]["profiles"][0, -1, :].astype(float)
        h = h - h.mean()
        ax.plot(np.arange(L), h, lw=0.7, label=rf"$\theta={th:.3f}$")
    ax.set_xlim(0, 400)
    ax.set_xlabel("$x$"); ax.set_ylabel("$h-\bar h$")
    ax.set_title(rf"(f) profiles, $L={L}$")
    ax.legend(fontsize=6)
    for a in axes.ravel()[:5]:
        a.axvline(THETA_C_RG, color="grey", ls="--", lw=0.9)
        a.axvline(THETA_C_LIT, color="grey", ls=":", lw=0.9)
    fig.tight_layout()
    save(fig, outdir, f"fig7_morphology_L{L}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main")
    ap.add_argument("--rule", default="bd")
    ap.add_argument("--skip-noise", action="store_true")
    args = ap.parse_args()

    outdir = os.path.join(ROOT, "figs", args.tag)
    os.makedirs(outdir, exist_ok=True)

    files = discover(args.tag, args.rule)
    cells = {k: load_result(v) for k, v in files.items()}
    print(f"{len(cells)} cells")

    adir = os.path.join(ROOT, "data", "analysis", f"{args.tag}_{args.rule}")
    rows = dict(np.load(os.path.join(adir, "rows.npz")))
    with open(os.path.join(adir, "summary.json")) as f:
        summary = json.load(f)

    if not args.skip_noise:
        fig_noise(outdir)
    fig_theta0(cells, outdir)
    Ls = sorted({L for (_, L) in cells})
    for L in (Ls[len(Ls) // 2], Ls[-1]):
        if (0.0, L) in cells:
            fig_widths(cells, outdir, L)
    fig_diagnostics(rows, outdir)
    fig_onset_drift(summary, outdir)
    fig_exponents(rows, outdir)
    for L in (Ls[len(Ls) // 2],):
        fig_morphology(rows, cells, outdir, L)
    print("all figures written")


if __name__ == "__main__":
    main()
