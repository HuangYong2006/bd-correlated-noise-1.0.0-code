"""
Morphology figure: what the interface actually looks like under the two
generators, at matched nominal theta.

Panels
 (a) representative interface profiles h(x) - <h>, stacked with offsets, for the
     corrected (spectral) generator;
 (b) the same for the FFGN generator;
 (c) the flat-terrace length distribution P(l) at L=256, both generators;
 (d) the structure factor S(k) at theta = 0 and 0.40, both generators.

Usage: python code/fig_morphology.py --L 256
"""

from __future__ import annotations

import argparse
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
    "axes.labelsize": 10, "axes.titlesize": 9.5, "legend.fontsize": 7,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.linewidth": 1.2,
})
TH_TH = 1.0 / 6.0
OUT = os.path.join(ROOT, "figs", "physics")

ap = argparse.ArgumentParser()
ap.add_argument("--L", type=int, default=256)
ap.add_argument("--rule", default="bd")
args = ap.parse_args()
L = args.L


def path(tag, theta):
    d = os.path.join(ROOT, "data", "raw", tag)
    tok = f"{theta:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
    p = os.path.join(d, f"{args.rule}_{tok}_L{L}.npz")
    return p if os.path.exists(p) else None


def profile(tag, theta, run=0):
    p = path(tag, theta)
    if p is None:
        return None
    z = np.load(p, allow_pickle=False)
    h = z["profiles"][run, -1, :].astype(np.float64)
    return h - h.mean()


def terrace_hist(tag, theta, max_len=25):
    p = path(tag, theta)
    if p is None:
        return None
    z = np.load(p, allow_pickle=False)
    H = z["profiles"][:, -1, :].astype(np.int64)
    counts = np.zeros(max_len + 1)
    tot = 0
    for i in range(min(48, H.shape[0])):
        d = H[i] - np.roll(H[i], -1)
        b = (d == 0).astype(np.int8)
        if b.all() or not b.any():
            continue
        idx = np.flatnonzero(b == 0)
        lengths = np.diff(np.append(idx, idx[0] + H.shape[1])) - 1
        lengths = lengths[lengths > 0]
        if lengths.size == 0:
            continue
        counts += np.bincount(np.minimum(lengths, max_len),
                              minlength=max_len + 1)
        tot += lengths.size
    return counts / max(tot, 1)


def spec(tag, theta):
    p = path(tag, theta)
    if p is None:
        return None
    z = np.load(p, allow_pickle=False)
    H = z["profiles"][:, -1, :].astype(np.float64)
    H = H - H.mean(axis=-1, keepdims=True)
    F = np.fft.rfft(H, axis=-1)
    S = (F.real ** 2 + F.imag ** 2) / H.shape[-1]
    return S[:, 1:H.shape[-1] // 2 + 1].mean(axis=0)


def main():
    os.makedirs(OUT, exist_ok=True)
    thetas = [0.0, TH_TH, 0.25, 0.40]
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(thetas)))
    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.0))

    for ax, tag, lab in ((axes[0, 0], "spec", "spectral (corrected)"),
                         (axes[0, 1], "main", "FFGN (published)")):
        off = 0.0
        for c, th in zip(cols, thetas):
            h = profile(tag, th, run=0)
            if h is None:
                continue
            off += 1.15 * h.std() + 6.0
            ax.plot(np.arange(L), h + off, lw=0.7, color=c,
                    label=rf"$\theta={th:.3f}$")
        ax.set_xlim(0, L)
        ax.set_xlabel("$x$")
        ax.set_ylabel(r"$h-\bar h$  (offset)")
        ax.set_title(lab)
        ax.legend(loc="upper right", fontsize=6.2)

    ax = axes[1, 0]
    for tag, mk, lab, col in (("spec", "o", "spectral", "C0"),
                              ("main", "s", "FFGN", "C3")):
        for th, ls in ((0.0, "-"), (0.40, "--")):
            v = terrace_hist(tag, th)
            if v is None:
                continue
            k = np.arange(v.size)
            ax.semilogy(k[1:], np.maximum(v[1:], 1e-6), mk + ls, ms=3,
                        color=col, alpha=1.0 if th == 0.0 else 0.6,
                        label=f"{lab}, $\\theta={th:.2f}$")
    ax.set_xlabel("terrace length $\\ell$  [lattice units]")
    ax.set_ylabel("$P(\\ell)$")
    ax.set_title(rf"(c) terrace-length distribution, $L={L}$")
    ax.legend(loc="upper right", fontsize=6.2)

    ax = axes[1, 1]
    for tag, lab, col in (("spec", "spectral", "C0"), ("main", "FFGN", "C3")):
        for th, ls in ((0.0, "-"), (0.40, "--")):
            S = spec(tag, th)
            if S is None:
                continue
            m = np.arange(1, S.size + 1)
            ax.loglog(m, S, ls, color=col, lw=1.0, alpha=1.0 if th == 0.0 else 0.7,
                      label=f"{lab}, $\\theta={th:.2f}$")
    ax.set_xlabel("$m$  ($k=2\\pi m/L$)")
    ax.set_ylabel("$S(k)$")
    ax.set_title(rf"(d) structure factor, $L={L}$")
    ax.legend(loc="lower left", fontsize=6.2)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_morphology.{ext}"), bbox_inches="tight")
    print("wrote figs/physics/fig_morphology.pdf/.png")

    # numeric summary
    print(f"\n{'tag':>6} {'theta':>8} {'sd(h)':>8} {'<l>':>8} {'alpha_s':>8}")
    for tag in ("spec", "main"):
        for th in thetas:
            h = profile(tag, th)
            v = terrace_hist(tag, th)
            S = spec(tag, th)
            if h is None:
                continue
            ell = float((np.arange(v.size) * v).sum())
            m = np.arange(1, S.size + 1)
            lo, hi = int(0.25 * S.size), int(0.5 * S.size)
            sl = np.polyfit(np.log(m[lo:hi]), np.log(S[lo:hi]), 1)[0]
            print(f"{tag:>6} {th:8.4f} {h.std():8.2f} {ell:8.3f} "
                  f"{-0.5*(sl+1):8.3f}")


if __name__ == "__main__":
    main()
