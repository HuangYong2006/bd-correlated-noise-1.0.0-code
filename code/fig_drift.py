r"""
Figure: does the deviation from KPZ scaling shrink with system size, and does
that conclusion survive when the noise filter has the EXACT fractional-Gaussian
correlation amplitude?

Three generators are compared at the two contested values of the parameter,
theta = 0.20 and theta = 0.25:

  FFGN          the published generator (Mandelbrot), whose binarised theta axis
                is offset by +0.19 and which never reaches the short-range regime
  powerlaw      spectral filtering with |f|^{-2 theta}; theta axis accurate to
                ~10%, under-correlates at large theta
  fgn           spectral filtering with the EXACT fGn spectral density; theta
                axis accurate to ~1%, exactly white at theta = 0

Panels
  (a) alpha_eff from consecutive size pairs versus L
  (b) the fixed-theta test: dH(theta, L) versus L.  A genuine crossover requires
      the deviation to GROW with L; L-independence means the observable carries
      no scale dependence.
  (c) alpha(theta) with bootstrap errors over runs
  (d) the split-half null floor versus L, the confound for any fixed-significance
      onset

Usage: python code/fig_drift.py
"""

from __future__ import annotations

import csv
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

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 9.5, "legend.fontsize": 6.6,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.linewidth": 1.25,
})
TH_TH = 1.0 / 6.0
OUT = os.path.join(ROOT, "figs", "physics")

# generator -> (tag, label, colour, linestyle)
GENS = [
    ("main", "FFGN", "C3", "--"),
    ("spec", "powerlaw", "C0", "-"),
    ("fgnspec", "exact fGn", "C2", "-."),
]
THETAS = [0.20, 0.25]
MARK = {0.20: "o", 0.25: "s"}


def have(tag):
    return os.path.exists(os.path.join(ROOT, "data", "analysis",
                                       f"{tag}_bd", "diagnostics.csv"))


def csv_cols(tag, rule="bd"):
    p = os.path.join(ROOT, "data", "analysis", f"{tag}_{rule}", "diagnostics.csv")
    with open(p, newline="") as f:
        r = list(csv.DictReader(f))
    out = {}
    for k in r[0]:
        try:
            out[k] = np.array([float(x[k]) for x in r])
        except (ValueError, KeyError):
            pass
    return out


def wsat_cells(tag, rule="bd"):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            z = np.load(f, allow_pickle=False)
            out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = \
                float(np.sqrt(z["w_series"].astype(np.float64).mean(axis=0))[-1])
    return out


def nearest(thetas, target):
    return min(thetas, key=lambda t: abs(t - target))


def alpha_pairs(tag, theta):
    W = wsat_cells(tag)
    Ls = sorted({L for (_, L) in W})
    th = nearest({t for (t, _) in W}, theta)
    xs, ys = [], []
    for a, b in zip(Ls[:-1], Ls[1:]):
        if (th, a) in W and (th, b) in W:
            xs.append(np.sqrt(a * b))
            ys.append(np.log(W[(th, b)] / W[(th, a)]) / np.log(b / a))
    return np.array(xs), np.array(ys)


def alpha_boot(tag, theta, n_boot=300, seed=0):
    d = os.path.join(ROOT, "data", "raw", tag)
    files = {}
    for f in sorted(glob.glob(os.path.join(d, "bd_*.npz"))):
        m = re.match(r"bd_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            if abs(round(sign * int(tok[3:]) / 10000.0, 4) - theta) < 1e-9:
                files[int(m.group("L"))] = f
    Ls = sorted(files)
    if len(Ls) < 3:
        return np.nan, np.nan
    cached = {L: np.load(files[L], allow_pickle=False)["w_series"].astype(np.float64)
              for L in Ls}
    W0 = [float(np.sqrt(cached[L].mean(axis=0))[-1]) for L in Ls]
    a0 = float(np.polyfit(np.log(Ls), np.log(W0), 1)[0])
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        v = [float(np.sqrt(cached[L][rng.integers(0, cached[L].shape[0],
                                                  cached[L].shape[0])]
                           .mean(axis=0))[-1]) for L in Ls]
        bs.append(np.polyfit(np.log(Ls), np.log(v), 1)[0])
    return a0, float(np.std(bs, ddof=1))


def main():
    os.makedirs(OUT, exist_ok=True)
    gens = [(t, l, c, ls) for t, l, c, ls in GENS if have(t)]
    if not gens:
        raise SystemExit("no analysis directories found; run code/analyze.py first")
    print("generators available:", [g[1] for g in gens])

    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.2))

    # ---- (a) alpha_eff from consecutive pairs -----------------------------
    ax = axes[0, 0]
    for tag, lab, col, ls in gens:
        for th in THETAS:
            xs, ys = alpha_pairs(tag, th)
            if xs.size:
                ax.semilogx(xs, ys, ls, marker=MARK[th], ms=3.5, color=col,
                            label=f"{lab}, $\\theta={th}$")
    ax.axhline(0.5, color="k", lw=1.0)
    ax.set_xlabel(r"$L$  (pair midpoint)")
    ax.set_ylabel(r"$\alpha_{\rm eff}$ from consecutive $L$ pairs")
    ax.set_title(r"(a) $\theta=0.20$ (circles) and $0.25$ (squares)")
    ax.legend(loc="upper right", ncol=2)

    # ---- (b) fixed-theta L-independence ----------------------------------
    ax = axes[0, 1]
    for tag, lab, col, ls in gens:
        d = csv_cols(tag)
        for th in THETAS:
            thn = nearest(set(d["theta"]), th)
            v, e, LL = [], [], []
            for L in sorted(set(d["L"].astype(int))):
                m = (d["L"].astype(int) == L) & (np.abs(d["theta"] - thn) < 1e-9)
                if m.sum():
                    v.append(d["dH"][m][0])
                    e.append(d["dH_err"][m][0] if "dH_err" in d else 0.0)
                    LL.append(L)
            if len(LL) >= 3:
                ax.errorbar(LL, v, yerr=e, fmt=MARK[th] + ls, ms=3.5, color=col,
                            label=f"{lab}, $\\theta={th}$", alpha=0.9)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.003)
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\Delta H(\theta,L)$ at fixed $\theta$")
    ax.set_title("(b) L-independent $\\Rightarrow$ no crossover resolved")
    ax.legend(loc="upper left", ncol=2)

    # ---- (c) alpha(theta) with bootstrap errors --------------------------
    ax = axes[1, 0]
    t = np.linspace(0.0, 0.45, 200)
    ax.plot(t, (1 + 4 * t) / (3 + 2 * t), "-", color="grey", lw=1.1,
            label="Medina one-loop LR")
    ax.plot(t, 1.43 * t + 0.26, ":", color="grey", lw=1.1,
            label="Squizzato--Canet NPRG")
    for tag, lab, col, ls in gens:
        ths = [0.0, 0.10, TH_TH, 0.20, 0.25, 0.30, 0.40]
        a, e, tt = [], [], []
        for th in ths:
            v, ev = alpha_boot(tag, th, 200)
            if np.isfinite(v):
                a.append(v); e.append(ev); tt.append(th)
        if a:
            ax.errorbar(tt, a, yerr=e, fmt=ls + "o", ms=4, color=col, label=lab)
    ax.axhline(0.5, color="k", lw=1.0)
    ax.axvline(TH_TH, color="grey", ls="--", lw=1)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\alpha$")
    ax.set_title(r"(c) $\alpha(\theta)$, bootstrap errors over runs")
    ax.legend(loc="upper left")

    # ---- (d) null floor --------------------------------------------------
    ax = axes[1, 1]
    for tag, lab, col, ls in gens:
        d = csv_cols(tag)
        Ls = sorted(set(d["L"].astype(int)))
        ax.errorbar(Ls, [d["Hnull"][d["L"].astype(int) == L][0] for L in Ls],
                    yerr=[d["Hnull_std"][d["L"].astype(int) == L][0] for L in Ls],
                    fmt=ls + "o", ms=3.5, color=col, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$L$")
    ax.set_ylabel("split-half null floor of $\\Delta H$")
    ax.set_title("(d) floor shrinks with $L$: confound for onset drifts")
    ax.legend(loc="upper right")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_drift.{ext}"), bbox_inches="tight")
    print("wrote figs/physics/fig_drift.pdf/.png")

    # ---- numeric summary -------------------------------------------------
    print(f"\n{'generator':>12} {'theta':>7} {'alpha':>9} {'err':>8} "
          f"{'alpha_eff pairs':>34}")
    for tag, lab, _c, _ls in gens:
        for th in THETAS:
            a, e = alpha_boot(tag, th, 300)
            xs, ys = alpha_pairs(tag, th)
            pr = ", ".join(f"{y:.3f}" for y in ys)
            print(f"{lab:>12} {th:7.2f} {a:9.4f} {e:8.4f}   {pr}")


if __name__ == "__main__":
    main()
