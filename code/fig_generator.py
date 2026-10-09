"""
Figure for the generator-bias result.

Panels:
 (a) normalised correlation C(tau)/C(0) of the FFGN noise (solid) against the
     exact fractional-Gaussian-noise form at the SAME nominal theta (dashed)
     and at the calibrated effective theta (dotted);
 (b) correlation amplitude: A_FFGN vs the exact A_fGn = 2 theta (theta+1/2),
     with the ratio on a twin axis;
 (c) the calibration theta_eff(theta_nom), showing that the whole FFGN family
     lies above the renormalisation-group threshold theta_th = 1/6;
 (d) the model-level consequence: the single structural order parameter
     p_flat = frac{Delta h = 0} versus theta, for the FFGN generator and for the
     theta-continuous spectral generator.

Cost is small (a few short FFGN sequences); the dynamics panels read existing
campaign cells.

Usage: python code/fig_generator.py
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
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from bd_core import ffgn_weights, ffgn_n_components, _u01

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.linewidth": 1.3,
})

TH_TH = 1.0 / 6.0
OUT = os.path.join(ROOT, "figs", "generator")


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


def fft_ac(x):
    n = x.size
    x = x - x.mean()
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    return ac / ac[0]


def exact_fgn(tau, theta_eff):
    H = theta_eff + 0.5
    t = tau.astype(np.float64)
    return 0.5 * ((t + 1) ** (2 * H) - 2 * t ** (2 * H) + np.abs(t - 1) ** (2 * H))


def theta_eff_from_corr(ac, taus):
    from scipy.optimize import minimize_scalar
    w = 1.0 / taus

    def cost(th):
        if th < 0 or th > 0.498:
            return 1e9
        return float(np.sum(w * (exact_fgn(taus, th) - ac[taus]) ** 2))

    r = minimize_scalar(cost, bounds=(0.0, 0.498), method="bounded",
                        options={"xatol": 1e-5})
    return float(r.x)


def pflat_from_cell(path):
    z = np.load(path, allow_pickle=False)
    h = z["profiles"][:, -1, :].astype(np.int64)
    d = h - np.roll(h, -1, axis=-1)
    return float((d == 0).mean())


def cells_by_theta(tag, rule="bd", L=256):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if not m or int(m.group("L")) != L:
            continue
        tok = m.group("tok")
        sign = -1.0 if tok[0] == "m" else 1.0
        try:
            out[round(sign * int(tok[3:]) / 10000.0, 4)] = pflat_from_cell(f)
        except Exception:
            pass
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    T, NREAL = 1 << 19, 3
    taus_fit = np.arange(1, 4097)
    tau = np.arange(1, 4097)

    nom = [0.0, 0.05, TH_TH, 0.25, 0.45]
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(nom)))

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.0))
    ax = axes[0, 0]
    rec = []
    for c, th in zip(cols, nom):
        N = ffgn_n_components(T)
        Wn, rn = ffgn_weights(th, N)
        if th == 0.0:
            Wn, rn = ffgn_weights(1e-9, N)
        ac = np.mean([fft_ac(ffgn_seq(T, Wn, rn, seed=300 + i))
                      for i in range(NREAL)], axis=0)
        te = theta_eff_from_corr(ac, taus_fit)
        rec.append((th, te))
        ax.loglog(tau, np.abs(ac[taus_fit]), color=c,
                  label=rf"FFGN $\theta_{{\rm nom}}={th:.3f}$")
        ax.loglog(tau, exact_fgn(tau, th), "--", color=c, lw=0.8, alpha=0.85)
        ax.loglog(tau, exact_fgn(tau, te), ":", color=c, lw=0.8, alpha=0.85)
    ax.set_xlabel(r"lag $\tau$  [sweeps]")
    ax.set_ylabel(r"$C(\tau)/C(0)$")
    ax.set_title("(a) FFGN vs exact fGn")
    ax.legend(loc="lower left", ncol=1)
    ax.plot([], [], "k--", lw=0.9, label="exact fGn, same $\\theta_{\\rm nom}$")
    ax.plot([], [], "k:", lw=0.9, label="exact fGn, $\\theta_{\\rm eff}$")
    ax.legend(loc="lower left", fontsize=6.2)
    rec = np.array(rec)

    # ---- (b) amplitudes ----
    try:
        amp = np.load(os.path.join(ROOT, "data", "generator_amplitude.npz"))["rows"]
    except Exception:
        amp = None
    ax = axes[0, 1]
    if amp is not None:
        m = amp[:, 0] <= 0.45
        tn = amp[m, 0]
        ax.semilogy(tn, np.maximum(amp[m, 1], 1e-3), "o-", ms=4, color="C3",
                    label=r"$A_{\rm FFGN}$ (measured)")
        ax.semilogy(tn, np.maximum(amp[m, 2], 1e-3), "s--", ms=4, color="C0",
                    label=r"$A_{\rm fGn}=2\theta(\theta+0.5)$")
        ax.semilogy(tn, np.maximum(amp[m, 3], 1e-3), "^:", ms=4, color="C2",
                    label=r"$A_{\rm spec}$ (spectral family)")
    ax.set_xlabel(r"nominal $\theta$")
    ax.set_ylabel(r"correlation amplitude $A$")
    ax.set_title(r"(b) amplitude: FFGN does not vanish at $\theta\to0$")
    ax.legend(loc="lower right")

    # ---- (c) calibration ----
    ax = axes[1, 0]
    try:
        pairs = np.load(os.path.join(ROOT, "data", "generator_calibration.npz"))["pairs"]
    except Exception:
        pairs = np.column_stack([rec[:, 0], rec[:, 1]])
    ax.plot([0, 0.45], [0, 0.45], "k--", lw=1, label="no distortion")
    ax.plot(pairs[:, 0], pairs[:, 1], "o-", ms=4, color="C3")
    ax.axhline(TH_TH, color="C0", ls="-", lw=1.2, label=r"$\theta_{\rm th}=1/6$")
    ax.axhspan(TH_TH, 0.5, color="C0", alpha=0.08, lw=0)
    ax.set_xlabel(r"nominal FFGN $\theta$")
    ax.set_ylabel(r"effective $\theta$")
    ax.set_ylim(0, 0.5)
    ax.set_title("(c) the FFGN family never reaches the SR-KPZ regime")
    ax.legend(loc="lower right")

    # ---- (d) model-level consequence ----
    ax = axes[1, 1]
    ff = cells_by_theta("main")
    sp = cells_by_theta("spec") or cells_by_theta("testspec")
    if ff:
        k = sorted(ff)
        ax.plot(k, [ff[t] for t in k], "o-", ms=4, color="C3",
                label="FFGN (published generator)")
    if sp:
        k = sorted(sp)
        ax.plot(k, [sp[t] for t in k], "s-", ms=4, color="C0",
                label="spectral ($\\theta$-continuous)")
    ax.axvline(TH_TH, color="grey", ls="--", lw=1)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$p_{\rm flat}={\rm frac}\{\Delta h=0\}$")
    ax.set_title(r"(d) $L=256$: spurious jump and wrong trend")
    ax.legend(loc="lower left")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_generator.{ext}"), bbox_inches="tight")
    print(f"wrote figs/generator/fig_generator.pdf/.png")
    print("theta_nom -> theta_eff:")
    for th, te in rec:
        print(f"   {th:7.4f} -> {te:7.4f}")


if __name__ == "__main__":
    main()
