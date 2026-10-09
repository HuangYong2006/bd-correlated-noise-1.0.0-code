r"""
Publication figures for the CTP manuscript.

FIGURE CONTRACT (established before any code, per the figure skill)
====================================================================
Journal/export contract: CTP is a SINGLE-COLUMN A4 journal with a text width of
178 mm (measured from a published CTP paper: 504 pt).  Every figure is therefore
exported at exactly 178 mm wide, PDF with editable TrueType text (fonttype 42),
plus 600 dpi TIFF and SVG.  Base font 7 pt, hard 5 pt floor for every glyph.

Palette: one restrained scheme -- neutral greys for reference material, a blue
family for the corrected/trusted results, a vermillion accent for the FFGN
generator that the paper argues is defective.  No twin axes, no filled colour
blocks, top/right spines removed, legend only where a direct label is impossible.

Archetype: quantitative grid, one hero panel per figure plus subordinate
evidence panels.

  FIG 1  "The FFGN generator's theta axis never intersects the short-range KPZ
          regime; the corrected filter is unbiased and exactly white at 0."
    (a) HERO  theta_eff(theta_nom) for both generators vs the identity, with the
              theta_th = 1/6 threshold band.  This is the whole argument.
    (b) amplitude ratio A_FFGN / A_fGn vs theta (why the axis is offset).
    (c) p_flat(theta) -- the consequence inside the model, sign inverted.

  FIG 2  "The measured roughness exponent is generator-dominated; the two
          corrected filters agree with each other."
    (a) HERO  alpha(theta) with bootstrap errors, three generators, against the
              Medina one-loop and Squizzato-Canet long-range branches.
    (b) alpha_eff from consecutive L pairs at theta = 0.25.

  FIG 3  "The fixed-theta scaling test detects a crossover only for theta >= 0.30."
    (a) HERO  d(dH)/dlnL vs theta, with significance marked.
    (b) dH(theta, L) vs L at four theta.

  FIG 4  "The interface morphology inverts between the two generators."
    (a) profiles at theta = 0.4 stacked; (b) terrace-length pmf; (c) S(k).

Usage: python code/figures_ctp.py
"""

from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

# Optional: the panel-alignment / collision audit helpers live in an external
# "nature-figure" skill directory.  It is NOT required to render the figures --
# if it is absent the audit gate is reported as skipped and the export proceeds.
# Point DSH_FIGURE_SKILL at your own copy, or just ignore the skip messages.
SKILL = os.environ.get("DSH_FIGURE_SKILL",
                       os.path.expanduser(r"~\.dsh\skills\nature-figure"))
if os.path.isdir(os.path.join(SKILL, "scripts")):
    sys.path.insert(0, os.path.join(SKILL, "scripts"))

# ----------------------------------------------------------------- export spec
MM = 1.0 / 25.4
# The saved PDF uses bbox="tight", so the page width is the ink extent plus
# 2*PAD_MM.  Sizing the canvas at (178 - 2*PAD_MM) makes the delivered page
# exactly 178 mm wide while still leaving a real margin: the first render had a
# 0.67 mm left margin, which put the "(a)" tag hard against the trim.
PAD_MM = 1.4
PAD_IN = PAD_MM * MM
W_FULL = 178 * MM - 2 * PAD_IN    # CTP single-column text width
W_HALF = 86 * MM
BASE_FS = 7.0

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "font.size": BASE_FS,
    "axes.labelsize": BASE_FS + 0.5,
    "axes.titlesize": BASE_FS + 0.5,
    "xtick.labelsize": BASE_FS - 0.5,
    "ytick.labelsize": BASE_FS - 0.5,
    "legend.fontsize": BASE_FS - 0.8,
    "axes.titlesize": BASE_FS + 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.6, "ytick.major.size": 2.6,
    "lines.linewidth": 1.1,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "legend.labelspacing": 0.40,
    "legend.borderaxespad": 0.3,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": PAD_IN,
})

# ------------------------------------------------------------------- palette
NEUTRAL = "#4d4d4d"
GREY_L = "#b3b3b3"
BLUE = "#1f6fb4"        # signal family: the corrected results
BLUE_L = "#7fb3d5"
VERM = "#d1495b"        # accent: the FFGN generator under critique
GREEN = "#3f8f6b"
BAND = "#d5e0ec"
OUT = os.path.join(ROOT, "figs", "ctp")
os.makedirs(OUT, exist_ok=True)

TH_TH = 1.0 / 6.0


# ------------------------------------------------------------------- data
def cells(tag, rule="bd"):
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, f"{rule}_*.npz"))):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if m:
            tok = m.group("tok")
            sign = -1.0 if tok[0] == "m" else 1.0
            z = np.load(f, allow_pickle=False)
            out[(round(sign * int(tok[3:]) / 10000.0, 4), int(m.group("L")))] = (
                float(np.sqrt(z["w_series"].astype(np.float64).mean(axis=0))[-1]),
                z["profiles"])
    return out


def pflat(tag, L=256):
    out = {}
    for (th, LL), (_, prof) in cells(tag).items():
        if LL != L:
            continue
        h = prof[:, -1, :].astype(np.int64)
        d = h - np.roll(h, -1, axis=-1)
        out[th] = float((d == 0).mean())
    return out


def cell_files(tag, theta, rule="bd"):
    """Map L -> path for one (tag, theta)."""
    d = os.path.join(ROOT, "data", "raw", tag)
    out = {}
    for f in glob.glob(os.path.join(d, f"{rule}_*.npz")):
        m = re.match(rf"{rule}_(?P<tok>[pm]0d\d+)_L(?P<L>\d+)\.npz$", os.path.basename(f))
        if not m:
            continue
        tok = m.group("tok")
        sign = -1.0 if tok[0] == "m" else 1.0
        if abs(round(sign * int(tok[3:]) / 10000.0, 4) - round(theta, 4)) < 1e-9:
            out[int(m.group("L"))] = f
    return out


def alpha_boot(tag, theta, n_boot=250, seed=0):
    """alpha from W_sat ~ L^alpha, bootstrap over runs.

    The per-cell w_series are loaded ONCE and reused across bootstrap draws;
    loading inside the loop costs thousands of compressed-archive reads.
    """
    files = cell_files(tag, theta)
    Ls = sorted(files)
    if len(Ls) < 3:
        return np.nan, np.nan
    W = {}
    for L in Ls:
        W[L] = np.load(files[L], allow_pickle=False)["w_series"].astype(np.float64)
    a0 = float(np.polyfit(np.log(Ls), np.log([np.sqrt(W[L].mean(axis=0))[-1]
                                             for L in Ls]), 1)[0])
    rng = np.random.default_rng(seed)
    bs = np.empty(n_boot)
    for k in range(n_boot):
        v = []
        for L in Ls:
            w = W[L]
            v.append(np.sqrt(w[rng.integers(0, w.shape[0], w.shape[0])]
                             .mean(axis=0))[-1])
        bs[k] = np.polyfit(np.log(Ls), np.log(v), 1)[0]
    return a0, float(bs.std(ddof=1))


def alpha_pairs(tag, theta):
    cs = cells(tag)
    Ls = sorted({L for (t, L) in cs if t == theta})
    xs, ys = [], []
    for a, b in zip(Ls[:-1], Ls[1:]):
        if (theta, a) in cs and (theta, b) in cs:
            xs.append(np.sqrt(a * b))
            ys.append(np.log(cs[(theta, b)][0] / cs[(theta, a)][0]) / np.log(b / a))
    return np.array(xs), np.array(ys)


def alpha_pairs_boot(tag, theta, n_boot=200, seed=0):
    """Pair-difference alpha_eff with a run-bootstrap error.

    Panel (b) of Fig. 2 shows the same three series as panel (a), so it must
    carry uncertainty too -- a bare three-point V reads as stronger evidence
    than it is.  The bootstrap resamples runs within each L, exactly as
    alpha_boot does, and propagates through both members of a pair.
    """
    files = cell_files(tag, theta)
    Ls = sorted(files)
    if len(Ls) < 2:
        return np.array([]), np.array([]), np.array([])
    W = {L: np.load(files[L], allow_pickle=False)["w_series"].astype(np.float64)
         for L in Ls}
    rng = np.random.default_rng(seed)
    xs, ys, es = [], [], []
    for a, b in zip(Ls[:-1], Ls[1:]):
        wa, wb = W[a], W[b]
        val = lambda w: np.sqrt(w.mean(axis=0))[-1]
        y = np.log(val(wb) / val(wa)) / np.log(b / a)
        ia = rng.integers(0, wa.shape[0], (n_boot, wa.shape[0]))
        ib = rng.integers(0, wb.shape[0], (n_boot, wb.shape[0]))
        ya = np.sqrt(wa[ia].mean(axis=1))[:, -1]
        yb = np.sqrt(wb[ib].mean(axis=1))[:, -1]
        bs = np.log(yb / ya) / np.log(b / a)
        xs.append(np.sqrt(a * b)); ys.append(y); es.append(float(bs.std(ddof=1)))
    return np.array(xs), np.array(ys), np.array(es)


def csv_cols(tag, rule="bd"):
    p = os.path.join(ROOT, "data", "analysis", f"{tag}_{rule}", "diagnostics.csv")
    with open(p, newline="") as f:
        r = list(__import__("csv").DictReader(f))
    out = {}
    for k in r[0]:
        try:
            out[k] = np.array([float(x[k]) for x in r])
        except (ValueError, KeyError):
            pass
    return out


def save(fig, name):
    """Panel-alignment gate, then export."""
    # An invisible full-canvas rectangle pins the "tight" bounding box to the
    # canvas.  Without it the page width follows the ink extent, which differed
    # by up to 7 mm between the four figures -- and since the manuscript scales
    # each figure to \\linewidth, that made the apparent font size differ
    # figure to figure.  With the box pinned, every page is exactly 178 mm.
    fig.add_artist(mpl.patches.Rectangle((0, 0), 1, 1,
                                         transform=fig.transFigure,
                                         fill=False, edgecolor="none",
                                         linewidth=0.0))
    fig.canvas.draw()
    tb = fig.get_tightbbox(fig.canvas.get_renderer())   # inches
    fw, fh = fig.get_size_inches()
    over = [f"{k} by {v*25.4:.2f} mm" for k, v in
            (("width", tb.width - fw), ("height", tb.height - fh)) if v > 1e-3]
    if over:
        print(f"    [canvas] ink exceeds the canvas in {', '.join(over)} "
              f"-- the page will be wider than 178 mm")
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(
            fig, json_out=os.path.join(OUT, f"{name}.alignment.json"),
            overlay_svg=os.path.join(OUT, f"{name}.alignment.svg"),
            tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    except Exception as e:
        print(f"    [alignment gate] {type(e).__name__}: {e}")
    fig.savefig(os.path.join(OUT, f"{name}.pdf"))
    fig.savefig(os.path.join(OUT, f"{name}.svg"))
    fig.savefig(os.path.join(OUT, f"{name}.tiff"), dpi=600)
    # preview raster, for visual review only -- not part of the submission set
    prev = os.environ.get("BD_FIG_PREVIEW", os.path.join(ROOT, ".figpreview"))
    os.makedirs(prev, exist_ok=True)
    fig.savefig(os.path.join(prev, f"{name}.png"), dpi=420)
    plt.close(fig)
    print(f"    wrote figs/ctp/{name}.pdf/.svg/.tiff  (+ preview raster)")


def panel_label(ax, s, dx=2.0, dy=4.0):
    """Bold (a)/(b)/(c) just above the top-left corner of the axes.

    The offset is in POINTS, not axes fractions.  An axes-fraction offset made
    the tag drift with the panel size: in Fig. 1 it sat 9 mm left of its own
    spine while its siblings sat 2 mm left, and in Fig. 4 the tag landed on top
    of the topmost y tick label.  A fixed 2 pt inset also keeps the tag clear of
    the tick-label column, which is what a tag above the spine collides with.
    """
    ax.annotate(s, xy=(0.0, 1.0), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=BASE_FS + 0.7, fontweight="bold",
                va="bottom", ha="left", annotation_clip=False)


def top_legend(ax, ncol=1, y=1.015, x=0.13, title=None):
    """Legend in a reserved strip ABOVE the axes.

    Placing the legend inside the axes is what made the first render fail the
    rendered collision audit: stroked data lines crossed the legend glyphs.
    A strip above the axes has no data in it, so nothing can cross it.  The
    strip is inset by ``x`` so the bold panel tag at x = -0.075 is separated
    from the first legend handle by a visible gap instead of reading as part
    of the same string.
    """
    ax.legend(loc="lower left", bbox_to_anchor=(x, y), ncol=ncol,
              borderaxespad=0.0, handletextpad=0.5, columnspacing=1.1,
              title=title,
              title_fontsize=BASE_FS - 0.6 if title else None)


# ============================================================ FIG 1
def fig1():
    print("  Fig 1: generator bias")
    fig = plt.figure(figsize=(W_FULL, 3.15))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0],
                          hspace=0.75, wspace=0.30,
                          left=0.07, right=0.99, top=0.86, bottom=0.11)

    # (a) HERO
    ax = fig.add_subplot(gs[0, :])
    cal = np.load(os.path.join(ROOT, "data", "generator_calibration_binarised.npz"))
    # The calibration arrays also contain the negative-theta control points.  They
    # lie outside this panel's xlim, and plotting them makes the polyline geometry
    # span the whole figure -- which the collision audit correctly flagged as
    # crossing the axis labels.  Plot only the range shown.
    ff = cal["ffgn"]; sp = cal["spectral"]
    ff = ff[(ff[:, 0] >= 0.0) & (ff[:, 0] <= 0.401)]
    sp = sp[(sp[:, 0] >= 0.0) & (sp[:, 0] <= 0.401)]
    x = np.linspace(0, 0.40, 10)
    band = ax.axhspan(TH_TH, 0.47, color=BAND, lw=0, zorder=0)
    ax.plot(x, x, "--", color=GREY_L, lw=0.8, zorder=1,
            label=r"unbiased: $\theta_{\rm eff}=\theta$")
    ax.plot(sp[:, 0], sp[:, 1], "s-", ms=3.2, color=BLUE, zorder=3,
            label="exact fGn filter (this work)")
    ax.plot(ff[:, 0], ff[:, 1], "o-", ms=3.4, color=VERM, zorder=4,
            label="FFGN (published work)")
    # the shaded region is the whole point of the panel, so it needs its own key
    ax.plot([], [], "s", ms=5.5, color=BAND,
            label=r"$\theta_{\rm eff}>\theta_{\rm th}$")
    ax.axhline(TH_TH, color=NEUTRAL, lw=0.8, ls="--", zorder=2)
    ax.text(0.398, TH_TH + 0.012, r"$\theta_{\rm th}=1/6$", fontsize=6,
            color=NEUTRAL, ha="right", va="bottom")
    # The vertical gap between the two generators at theta_nom = 0.  It sits in
    # the one region of this panel that no polyline enters (the wedge bounded
    # below by the unbiased line and above by the FFGN curve).
    ax.annotate("", xy=(0.009, 0.1888), xytext=(0.009, 0.0),
                arrowprops=dict(arrowstyle="<->", color=VERM, lw=0.7))
    ax.text(0.018, 0.098, r"$+0.19$", fontsize=6.5, color=VERM,
            ha="left", va="center")
    ax.set_xlim(-0.008, 0.408); ax.set_ylim(0, 0.47)
    ax.set_xlabel(r"nominal $\theta$")
    ax.set_ylabel(r"effective $\theta$", labelpad=2.0)
    ax.tick_params(pad=1.8)
    panel_label(ax, "(a)")
    top_legend(ax, ncol=2, y=1.02)

    # (b) amplitude ratio
    ax = fig.add_subplot(gs[1, 0])
    amp = np.load(os.path.join(ROOT, "data", "generator_amplitude.npz"))["rows"]
    # the exact fGn amplitude is 2*theta*(theta+1/2), which vanishes at theta=0,
    # so the ratio is genuinely infinite there and theta=0 is left out
    m = (amp[:, 0] > 0.0) & (amp[:, 0] <= 0.401)
    r = amp[m, 1] / amp[m, 2]
    ax.semilogy(amp[m, 0], r, "o-", ms=3.0, color=VERM)
    ax.axhline(1.0, color=NEUTRAL, lw=0.8, ls="--")
    # the label used to sit on the dashed line and read as part of the dashes;
    # it is also kept well clear of the bottom spine, which the audit counts as
    # a stroked path crossing the text box
    ax.text(0.055, 0.60, r"exact $A_{\rm fGn}$", fontsize=6, color=NEUTRAL,
            ha="left", va="center")
    # every ratio is > 1, so there is no decade below the data; the default log
    # minor ticks packed a solid black comb against the spine
    ax.set_ylim(0.30, None)
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(-0.008, 0.408)
    ax.set_xlabel(r"nominal $\theta$")
    ax.set_ylabel(r"$A_{\rm FFGN}/A_{\rm fGn}$")
    panel_label(ax, "(b)")

    # (c) p_flat
    ax = fig.add_subplot(gs[1, 1])
    for tag, lab, col, mk, z in (("spec", "exact fGn filter", BLUE, "s", 3),
                                 ("main", "FFGN", VERM, "o", 4)):
        d = pflat(tag)
        k = sorted(t for t in d if 0.0 <= t <= 0.401)
        # zorder keeps the FFGN circle on top: in the first render the blue
        # square sat exactly on the red point at theta = 0 and hid it, leaving
        # a stray red stub with no visible endpoint
        ax.plot(k, [d[t] for t in k], mk + "-", ms=3.0, color=col, label=lab,
                zorder=z)
    ax.set_xlim(-0.008, 0.408)
    ax.set_xlabel(r"nominal $\theta$")
    ax.set_ylabel(r"$p_{\rm flat}$")
    panel_label(ax, "(c)")
    top_legend(ax, ncol=1, y=1.02)
    save(fig, "fig1_generator")


# ============================================================ FIG 2
def fig2():
    print("  Fig 2: exponents")
    fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 2.95),
                             gridspec_kw=dict(wspace=0.26, left=0.07,
                                              right=0.99, top=0.80, bottom=0.14))

    ths = [0.0, 0.10, TH_TH, 0.20, 0.25, 0.30, 0.40]
    t = np.linspace(0, 0.45, 200)
    ax = axes[0]
    ax.plot(t, (1 + 4 * t) / (3 + 2 * t), "-", color=GREY_L, lw=1.0,
            label="Medina (one-loop)")
    ax.plot(t, 1.43 * t + 0.26, ":", color=GREY_L, lw=1.1,
            label="Squizzato--Canet")
    ax.axhline(0.5, color=NEUTRAL, lw=0.8, label=r"$\alpha=1/2$")
    # the guide line stops BELOW the label's band: if it ran past the text the
    # audit counts it as a stroked path crossing a text box
    ax.axvline(TH_TH, color=GREY_L, lw=0.7, ls="--", zorder=0, ymax=0.80)
    ax.text(TH_TH + 0.014, 0.885, r"$\theta_{\rm th}=1/6$", fontsize=6,
            color=NEUTRAL, ha="left", va="center")
    for tag, lab, col, mk, filled in (("main", "FFGN", VERM, "o", True),
                                      ("spec", "power law", BLUE_L, "s", False),
                                      ("fgnspec", "exact fGn", BLUE, "^", True)):
        a, e, tt = [], [], []
        for th in ths:
            v, ev = alpha_boot(tag, round(th, 4), 250)
            if np.isfinite(v):
                a.append(v); e.append(ev); tt.append(th)
        # The two corrected filters agree to within 0.008, so coincident markers
        # are the expected outcome.  They are drawn as marker + separate bars
        # (which also keeps error-bar caps out of the legend handle) and the
        # power-law series is left OPEN so the two never hide each other.
        kw = dict(ms=3.6 if not filled else 3.2)
        if not filled:
            kw.update(markerfacecolor="white", markeredgewidth=0.9)
        ax.plot(tt, a, mk + "-", color=col, label=lab, zorder=3, **kw)
        ax.errorbar(tt, a, yerr=e, fmt="none", ecolor=col, capsize=1.3,
                    elinewidth=0.7, zorder=2)
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\alpha$")
    panel_label(ax, "(a)")
    top_legend(ax, ncol=2, y=1.02)

    ax = axes[1]
    for tag, lab, col, mk in (("main", "FFGN", VERM, "o"),
                              ("spec", "power law", BLUE_L, "s"),
                              ("fgnspec", "exact fGn", BLUE, "^")):
        xs, ys, es = alpha_pairs_boot(tag, 0.25, 200)
        if xs.size:
            ax.semilogx(xs, ys, mk + "-", ms=3.2, color=col, label=lab)
            ax.errorbar(xs, ys, yerr=es, fmt="none", ecolor=col, capsize=1.3,
                        elinewidth=0.7)
    ax.axhline(0.5, color=NEUTRAL, lw=0.8, label=r"$\alpha=1/2$")
    ax.set_ylim(0.46, None)
    # plain tick labels: the default 6x10^1 / 2x10^2 forms put ~4 pt superscript
    # exponents at the right-hand end of the axis
    ax.set_xticks([50, 100, 200, 400])
    ax.set_xticklabels(["50", "100", "200", "400"])
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    # the pair midpoints run from sqrt(32*64)=45 to sqrt(256*512)=362; the first
    # render set the limit AT 45 and bisected that marker on the left spine
    ax.set_xlim(38, 450)
    ax.set_xlabel(r"$L$ (pair midpoint, lattice sites)")
    ax.set_ylabel(r"$\alpha_{\rm eff}$")
    panel_label(ax, "(b)")
    top_legend(ax, ncol=2, y=1.02)
    save(fig, "fig2_exponents")


# ============================================================ FIG 3
def fig3():
    print("  Fig 3: fixed-theta scaling test")
    d_sp = csv_cols("spec"); d_fg = csv_cols("fgnspec")
    fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 2.95),
                             gridspec_kw=dict(wspace=0.26, left=0.085,
                                              right=0.99, top=0.80, bottom=0.12))

    ax = axes[0]
    for d, lab, col, mk, ls, filled in (
            (d_sp, "power law", BLUE_L, "s", "--", False),
            (d_fg, "exact fGn", BLUE, "^", "-", True)):
        xs, ys, sig = [], [], []
        for th in sorted(set(np.round(d["theta"], 4))):
            m = np.abs(d["theta"] - th) < 1e-9
            if not m.any() or th < 0:
                continue
            L, v, e = d["L"][m], d["dH"][m], d["dH_err"][m]
            o = np.argsort(L); L, v, e = L[o], v[o], e[o]
            if len(L) < 3:
                continue
            sl = np.polyfit(np.log(L), v, 1)[0]
            xs.append(th); ys.append(sl)
            sig.append(abs(sl) * np.log(L[-1] / L[0]) > 2 * np.mean(e))
        # plotted in units of 1e-3: the raw values span +-0.007, and tick labels
        # like "-0.0060" are wide enough to push the y-label off the canvas
        kw = dict(ms=3.4 if filled else 3.6)
        if not filled:
            kw.update(markerfacecolor="white", markeredgewidth=0.9)
        ax.plot(xs, np.array(ys) * 1e3, mk + ls, color=col, label=lab,
                zorder=3, **kw)
        # Significance is marked with an OPEN RING around the point.  The first
        # render stamped a filled star concentric with the marker, which hid the
        # square/triangle shape at exactly the points the panel is about.
        fx = [x for x, s in zip(xs, sig) if s]
        fy = [y * 1e3 for y, s in zip(ys, sig) if s]
        if fx:
            # Two distinct ring radii, because at theta >= 0.30 BOTH filters are
            # significant.  The radii differ by 2.1 pt while the worst measured
            # centre offset between the two series is 1.0 pt, so the rings stay
            # nested (1.1 pt of white between them) instead of merging.
            ax.plot(fx, fy, "o", ms=10.6 if lab.startswith("power") else 6.4,
                    markerfacecolor="none",
                    markeredgecolor=col, markeredgewidth=0.8, zorder=4,
                    label=r"$2\sigma$" if lab.startswith("power")
                    else None)
    ax.axhline(0, color=NEUTRAL, lw=0.8, zorder=1)
    # The threshold guide is keyed in the legend strip rather than annotated
    # inside the axes: the panel is dense enough that any in-axes position for
    # the label was crossed by a data curve in one render or another.
    ax.axvline(TH_TH, color=GREY_L, lw=0.7, ls="--", zorder=0, ymax=1.0,
               label=r"$\theta_{\rm th}=1/6$")
    lim = max(abs(np.array(ax.get_ylim())).max(), 1e-3) * 1.30
    ax.set_ylim(-lim, lim)          # headroom for the largest excursions
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$10^{3}\,\mathrm{d}\,\Delta H/\mathrm{d}\ln L$")
    panel_label(ax, "(a)")
    top_legend(ax, ncol=2, y=1.02)

    ax = axes[1]
    cmap = plt.cm.Blues(np.linspace(0.40, 0.95, 4))
    Lall = []
    for c, th in zip(cmap, [0.10, 0.20, 0.25, 0.30]):
        m = np.abs(d_fg["theta"] - round(th, 4)) < 1e-9
        L, v, e = d_fg["L"][m], d_fg["dH"][m], d_fg["dH_err"][m]
        o = np.argsort(L)
        Lall.append(L[o])
        ax.plot(L[o], v[o], "o-", ms=3.0, color=c, label=rf"$\theta={th}$")
        ax.errorbar(L[o], v[o], yerr=e[o], fmt="none", ecolor=c, capsize=1.2,
                    elinewidth=0.7)
    ax.axhline(0, color=NEUTRAL, lw=0.8)
    ax.set_xscale("log")
    Ls_all = np.unique(np.concatenate(Lall))
    ax.set_xticks(Ls_all)
    ax.set_xticklabels([f"{int(v)}" for v in Ls_all])
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(Ls_all.min() * 0.82, Ls_all.max() * 1.20)
    ax.set_xlabel(r"$L$ (lattice sites)")
    ax.set_ylabel(r"$\Delta H(\theta,L)$")
    panel_label(ax, "(b)")
    top_legend(ax, ncol=2, y=1.02, title="exact fGn filter")
    save(fig, "fig3_fixed_theta")


# ============================================================ FIG 4
def fig4():
    print("  Fig 4: morphology")
    L = 256
    cs_sp = cells("spec"); cs_ff = cells("main")
    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 2.75),
                             gridspec_kw=dict(wspace=0.30, left=0.05,
                                              right=0.99, top=0.80, bottom=0.155))

    ax = axes[0]
    # The first render of this panel was BLANK.  The offset accumulator ran
    # negative while the y-limits were computed from a sign-flipped copy of it,
    # so the window ended up inverted (bottom above top) and every profile fell
    # outside it.  Offsets are now accumulated with the sign that moves each
    # profile DOWN, and the limits are taken from the plotted arrays themselves.
    gap = 10.0
    y = 0.0
    lo, hi = np.inf, -np.inf
    for tag, lab, col in (("spec", "exact fGn filter", BLUE),
                          ("main", "FFGN", VERM)):
        prof = (cs_sp if tag == "spec" else cs_ff)[(0.40, L)][1]
        h = prof[0, -1, :].astype(float); h -= h.mean()
        ax.plot(np.arange(L), h + y, lw=0.5, color=col, label=lab, zorder=3)
        ax.axhline(y, color=GREY_L, lw=0.5, ls=":", zorder=1)
        lo = min(lo, float((h + y).min())); hi = max(hi, float((h + y).max()))
        y -= (np.ptp(h) + gap)
    # one h = 10 scale bar, placed ABOVE the topmost trace so no profile can
    # wander into it (the first attempt put it inside the trace's own band and
    # the trace ran straight through the "10")
    xb = L - 20
    yb = hi + 12.0
    ax.plot([xb, xb], [yb, yb - 10.0], color=NEUTRAL, lw=0.9, zorder=4)
    ax.text(xb - 5, yb - 5.0, "10", fontsize=6, color=NEUTRAL, ha="right",
            va="center", zorder=4)
    ax.set_xlim(0, L); ax.set_ylim(lo - 5.0, hi + 16.0)
    ax.set_xlabel(r"$x$ (lattice sites)")
    ax.set_ylabel(r"$h-\bar h$  (offset)")
    ax.set_yticks([])
    panel_label(ax, "(a)")
    top_legend(ax, ncol=1, y=1.02)

    ax = axes[1]
    for tag, lab, col, mk, mfc in (("spec", "exact fGn filter", BLUE, "o", BLUE),
                                   ("main", "FFGN", VERM, "s", "white")):
        prof = (cs_sp if tag == "spec" else cs_ff)[(0.40, L)][1]
        H = prof[:, -1, :].astype(np.int64)
        cnt = np.zeros(26); tot = 0
        for i in range(min(48, H.shape[0])):
            d = H[i] - np.roll(H[i], -1)
            b = (d == 0).astype(np.int8)
            if b.all() or not b.any():
                continue
            idx = np.flatnonzero(b == 0)
            ln = np.diff(np.append(idx, idx[0] + H.shape[1])) - 1
            ln = ln[ln > 0]
            if ln.size:
                cnt += np.bincount(np.minimum(ln, 25), minlength=26); tot += ln.size
        if tot:
            v = cnt / tot
            # Empty bins are left out rather than floored at an invented 1e-5:
            # the previous floor printed a flat plateau that looked like measured
            # structure.  The two series are also drawn with different shapes and
            # one open marker, because they coincide on that plateau and the
            # filled red circles used to hide the blue ones completely.
            yy = v[1:].astype(float)
            yy[yy <= 0] = np.nan
            ax.semilogy(np.arange(1, 26), yy, mk + "-", ms=2.8, color=col,
                        markerfacecolor=mfc, markeredgewidth=0.6,
                        label=lab, zorder=3 if mfc == "white" else 2)
    ax.set_xlabel(r"terrace length $\ell$"); ax.set_ylabel(r"$P(\ell)$")
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    panel_label(ax, "(b)")
    top_legend(ax, ncol=1, y=1.02)

    ax = axes[2]
    for tag, lab, col in (("spec", "exact fGn filter", BLUE),
                          ("main", "FFGN", VERM)):
        prof = (cs_sp if tag == "spec" else cs_ff)[(0.40, L)][1].astype(np.float64)
        H = prof[:, -1, :]; H = H - H.mean(axis=-1, keepdims=True)
        F = np.fft.rfft(H, axis=-1)
        S = (F.real ** 2 + F.imag ** 2).mean(axis=0) / H.shape[-1]
        S = S[1:H.shape[-1] // 2 + 1]
        ax.loglog(np.arange(1, S.size + 1), S, "-", lw=0.9, color=col, label=lab)
    ax.set_xlabel(r"$m$ ($k=2\pi m/L$)"); ax.set_ylabel(r"$S(k)$")
    panel_label(ax, "(c)")
    top_legend(ax, ncol=1, y=1.02)
    save(fig, "fig4_morphology")


if __name__ == "__main__":
    print(f"export width {W_FULL:.3f} in = {W_FULL*25.4:.0f} mm, base font {BASE_FS} pt")
    fig1(); fig2(); fig3(); fig4()
    print("done")
