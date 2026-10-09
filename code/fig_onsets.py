"""
Headline figure: the two calibrated diagnostics and the drift of the
scale-dependent onset, on the corrected (theta-continuous) axis.

The onsets are recomputed here directly from the stored per-cell table so that
bootstrap error bars are available for every criterion.

Usage: python code/fig_onsets.py [--tag spec] [--rule bd]
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from diagnostics import onset_with_error, fit_onset_scaling

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 7,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.linewidth": 1.3,
})
TH_TH = 1.0 / 6.0
OUT = os.path.join(ROOT, "figs", "physics")

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="spec")
ap.add_argument("--rule", default="bd")
args = ap.parse_args()

import csv

_adir = os.path.join(ROOT, "data", "analysis", f"{args.tag}_{args.rule}")
_csv = os.path.join(_adir, "diagnostics.csv")
with open(_csv, newline="") as f:
    _r = list(csv.DictReader(f))
rows = {}
for _k in _r[0]:
    try:
        rows[_k] = np.array([float(x[_k]) for x in _r], dtype=np.float64)
    except (ValueError, KeyError):
        pass
th = rows["theta"]
Ls = sorted(set(rows["L"].astype(int)))
print(f"{args.tag}: Ls = {Ls}  ({len(_r)} cells)")

fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.2))

# ---- (a)-(c): observables vs theta with null bands ----
for ax, L in zip(axes.ravel()[:3], Ls):
    m = rows["L"].astype(int) == L
    t = th[m]
    o = np.argsort(t)
    t = t[o]
    dH = rows["dH"][m][o]
    eH = rows["dH_err"][m][o]
    hn = rows["Hnull"][m][o]
    hn_sd = rows["Hnull_std"][m][o]
    phi = rows["dPhi"][m][o]

    ax.errorbar(t, dH, yerr=eH, fmt="o-", ms=3.5, color="C0",
                label=r"Hellinger $\Delta H$ (increment pmf)")
    ax.fill_between(t, -2 * hn_sd, 2 * hn_sd, color="grey", alpha=0.45, lw=0)
    ax.axhline(0, color="k", lw=0.7)
    ax.axvline(TH_TH, color="C3", ls="--", lw=1.1)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_title(rf"$L={L}$")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\Delta H$")
    ax2 = ax.twinx()
    ax2.plot(t, phi, "s--", ms=3, color="C1", lw=1.0,
             label=r"iForest $\Phi$")
    ax2.set_ylabel(r"iForest anomaly fraction $\Phi$", color="C1")
    ax2.tick_params(axis="y", colors="C1")
    ax2.grid(False)
    if L == Ls[0]:
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=6.5)
axes.ravel()[0].plot([], [], color="C3", ls="--", label=r"$\theta_{\rm th}=1/6$")
axes.ravel()[0].legend(loc="upper left", fontsize=6.2)
axes.ravel()[0].set_title(rf"(a) $L={Ls[0]}$", loc="left")

# ---- (d): onset drift ----
ax = axes.ravel()[3]
criteria = [
    ("Hellinger", "zH", [3.0, 5.0], ["o", "s"]),
    ("iForest", "dPhi", [0.10], ["^"]),
]
Ls_arr = np.array(Ls, dtype=float)
rec = []
for lab, key, thrs, mk in criteria:
    for zt, m1 in zip(thrs, mk):
        ts, es = [], []
        for L in Ls:
            m = rows["L"].astype(int) == L
            t = rows["theta"][m]
            v = rows[key][m].astype(np.float64)
            e = rows.get(key + "_err", np.full_like(v, np.nan))[m] \
                if (key + "_err") in rows else np.full_like(v, np.nan)
            e = np.where(np.isfinite(e), e, np.nanmedian(e) if np.isfinite(e).any() else 1.0)
            ob = onset_with_error(t, v, e, zt, True, 400, 3)
            ts.append(ob["theta_star"])
            es.append(ob["theta_star_err"])
        ts = np.array(ts, dtype=float)
        es = np.array(es, dtype=float)
        ok = np.isfinite(ts)
        if ok.sum() < 2:
            continue
        ax.errorbar(Ls_arr[ok], ts[ok], yerr=es[ok], fmt=m1 + "-", ms=5,
                    label=f"{lab}, crit {zt}")
        rec.append((lab, zt, ts, es))
        if ok.sum() >= 3:
            f0 = fit_onset_scaling(Ls_arr[ok], ts[ok],
                                   es[ok] if np.isfinite(es[ok]).all() else None,
                                   theta_c_free=False)
            if np.isfinite(f0["nu"]):
                Lg = np.geomspace(Ls_arr[ok].min(), Ls_arr[ok].max() * 4, 60)
                ax.plot(Lg, f0["A"] * Lg ** (-1.0 / f0["nu"]), ":", lw=1,
                        color=ax.lines[-1].get_color())
ax.axhline(TH_TH, color="C3", ls="--", lw=1.2, label=r"$\theta_{\rm th}=1/6$")
ax.axhline(0.0, color="k", lw=1)
ax.set_xscale("log")
ax.set_xlabel(r"$L$")
ax.set_ylabel(r"onset $\theta^*(L)$")
ax.set_title("(d) onset drift (dotted: forced-zero FSS fit)")
ax.legend(loc="upper right", fontsize=6.3)
axes.ravel()[3].set_title("(d) onset drift (dotted: forced-zero FSS fit)",
                          loc="left")

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"fig_onsets.{ext}"), bbox_inches="tight")
print("wrote figs/physics/fig_onsets.pdf/.png\n")

print(f"{'diagnostic':>12} {'crit':>6} " + " ".join(f"{'L='+str(L):>18}" for L in Ls))
for lab, zt, ts, es in rec:
    cells = " ".join(f"{t:8.4f}+-{e:.4f}" if np.isfinite(t) else f"{'--':>16}"
                     for t, e in zip(ts, es))
    print(f"{lab:>12} {zt:6.2f} {cells}")
