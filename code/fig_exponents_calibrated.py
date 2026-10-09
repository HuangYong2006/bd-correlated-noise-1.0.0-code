"""
Key physics figure: does the measured roughness exponent depend on the noise
GENERATOR, and is any difference explained by the generator's theta-axis offset?

Crucially, alpha is refitted here over a COMMON system-size window for both
tracks, because alpha is extracted from W_sat ~ L^alpha and finite-size
corrections make the fitted value depend on the L range used.

Panels
 (a) alpha vs NOMINAL theta for both generators
 (b) alpha vs EFFECTIVE (binarised) theta for both generators, with the
     renormalisation-group threshold marked
 (c) the binarised theta_eff(theta_nom) calibration of both generators
 (d) p_flat(theta) for both generators

Usage: python code/fig_exponents_calibrated.py
"""

from __future__ import annotations

import glob
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
from bd_core import load_result

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.linewidth": 1.3,
})
TH_TH = 1.0 / 6.0
COMMON_L = [64, 128, 256]
OUT = os.path.join(ROOT, "figs", "physics")


def cells(tag, rule="bd"):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = f
    return out


def wsat(path):
    try:
        z = np.load(path, allow_pickle=False)
        return float(np.sqrt(z["w_series"].astype(np.float64).mean(axis=0))[-1])
    except Exception:
        return np.nan


def alpha_beta(tag, Ls):
    cs = cells(tag)
    thetas = sorted({t for (t, L) in cs})
    out = {}
    for th in thetas:
        W, LL = [], []
        for L in Ls:
            if (th, L) in cs:
                w = wsat(cs[(th, L)])
                if np.isfinite(w):
                    W.append(w); LL.append(L)
        if len(LL) >= 3:
            p = np.polyfit(np.log(LL), np.log(W), 1)
            res = np.log(W) - np.polyval(p, np.log(LL))
            sxx = ((np.log(LL) - np.log(LL).mean()) ** 2).sum()
            err = float(np.sqrt((res ** 2).sum() / max(1, len(LL) - 2) / max(sxx, 1e-12)))
            out[th] = (float(p[0]), err, len(LL))
    return out


def pflat(tag, L=256):
    cs = cells(tag)
    out = {}
    for (th, LL), f in cs.items():
        if LL != L:
            continue
        try:
            z = np.load(f, allow_pickle=False)
            h = z["profiles"][:, -1, :].astype(np.int64)
            d = h - np.roll(h, -1, axis=-1)
            out[th] = float((d == 0).mean())
        except Exception:
            pass
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    ab_ff = alpha_beta("main", COMMON_L)
    ab_sp = alpha_beta("spec", COMMON_L)
    try:
        cal = np.load(os.path.join(ROOT, "data",
                                   "generator_calibration_binarised.npz"))
        ff_cal, sp_cal = cal["ffgn"], cal["spectral"]
    except Exception:
        ff_cal = sp_cal = None

    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.2))

    def series(ab):
        th = np.array(sorted(ab))
        a = np.array([ab[t][0] for t in th])
        e = np.array([ab[t][1] for t in th])
        return th, a, e

    ax = axes[0, 0]
    tf, af, ef = series(ab_ff)
    ts, as_, es = series(ab_sp)
    ax.errorbar(tf, af, yerr=ef, fmt="o-", ms=4, color="C3",
                label="FFGN (published generator)")
    ax.errorbar(ts, as_, yerr=es, fmt="s-", ms=4, color="C0",
                label=r"spectral ($\theta$-continuous)")
    ax.axhline(0.5, color="k", lw=0.9)
    ax.axvline(TH_TH, color="grey", ls="--", lw=1)
    ax.set_xlabel(r"nominal $\theta$"); ax.set_ylabel(r"$\alpha$")
    ax.set_title(rf"(a) $\alpha$ vs nominal $\theta$  ($L={COMMON_L}$)")
    ax.legend(loc="upper left")

    ax = axes[0, 1]
    if ff_cal is not None:
        etf = np.interp(tf, ff_cal[:, 0], ff_cal[:, 1])
        ets = np.interp(ts, sp_cal[:, 0], sp_cal[:, 1])
        ax.errorbar(etf, af, yerr=ef, fmt="o", ms=4, color="C3", label="FFGN")
        ax.errorbar(ets, as_, yerr=es, fmt="s", ms=4, color="C0",
                    label="spectral")
    ax.axhline(0.5, color="k", lw=0.9)
    ax.axvline(TH_TH, color="grey", ls="--", lw=1)
    ax.set_xlabel(r"effective $\theta$ (binarised, calibrated)")
    ax.set_ylabel(r"$\alpha$")
    ax.set_title(r"(b) $\alpha$ vs effective $\theta$")
    ax.legend(loc="upper left")

    ax = axes[1, 0]
    if ff_cal is not None:
        ax.plot(ff_cal[:, 0], ff_cal[:, 1], "o-", ms=4, color="C3",
                label="FFGN")
        ax.plot(sp_cal[:, 0], sp_cal[:, 1], "s-", ms=4, color="C0",
                label="spectral")
        ax.plot([0, 0.45], [0, 0.45], "k--", lw=1, label="no offset")
    ax.axhline(TH_TH, color="C2", lw=1.4)
    ax.axhspan(TH_TH, 0.5, color="C2", alpha=0.08, lw=0)
    ax.set_xlabel(r"nominal $\theta$")
    ax.set_ylabel(r"effective $\theta$")
    ax.set_ylim(0, 0.48)
    ax.set_title(r"(c) $\theta$ axis: FFGN never enters $\theta_{\rm eff}<1/6$")
    ax.legend(loc="lower right")

    ax = axes[1, 1]
    pf_ff = pflat("main"); pf_sp = pflat("spec")
    if pf_ff:
        k = sorted(pf_ff)
        ax.plot(k, [pf_ff[t] for t in k], "o-", ms=4, color="C3", label="FFGN")
    if pf_sp:
        k = sorted(pf_sp)
        ax.plot(k, [pf_sp[t] for t in k], "s-", ms=4, color="C0",
                label="spectral")
    ax.axvline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"nominal $\theta$")
    ax.set_ylabel(r"$p_{\rm flat}$")
    ax.set_title(r"(d) $L=256$: FFGN jumps and inverts the trend")
    ax.legend(loc="lower left")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_exponents_calibrated.{ext}"),
                    bbox_inches="tight")
    print("wrote figs/physics/fig_exponents_calibrated.pdf/.png")
    print(f"\ncommon L window used for alpha: {COMMON_L}")
    print(f"{'theta':>9} {'FFGN alpha':>12} {'spectral alpha':>15} "
          f"{'FFGN th_eff':>12} {'spec th_eff':>12}")
    for th in sorted(set(list(ab_ff) + list(ab_sp))):
        a1 = f"{ab_ff[th][0]:.4f}" if th in ab_ff else "  --  "
        a2 = f"{ab_sp[th][0]:.4f}" if th in ab_sp else "  --  "
        e1 = f"{np.interp(th, ff_cal[:,0], ff_cal[:,1]):.4f}" if ff_cal is not None else "--"
        e2 = f"{np.interp(th, sp_cal[:,0], sp_cal[:,1]):.4f}" if sp_cal is not None else "--"
        print(f"{th:9.4f} {a1:>12} {a2:>15} {e1:>12} {e2:>12}")


if __name__ == "__main__":
    main()
