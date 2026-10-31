"""Revised PRE figures. Read only the new summary and its declared raw files."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
from matplotlib.ticker import ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
STYLE = {"shared": ("Shared AR", "#d1495b", "o"),
         "powerlaw": ("Power law", "#7fb3d5", "s"),
         "fgn": ("fGn spectrum", "#1f6fb4", "^")}
NEUTRAL = "#4d4d4d"


def configure_style():
    mpl.rcParams.update({
        "font.family": "Times New Roman",
        "font.serif": ["Times New Roman"],
        "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 7.5,
        "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 6.2,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
        "mathtext.sf": "Times New Roman",
        "mathtext.tt": "Times New Roman",
        "mathtext.cal": "Times New Roman",
        "mathtext.fallback": None, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": 0.7,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.major.size": 2.6, "ytick.major.size": 2.6,
        "lines.linewidth": 1.1, "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "figure.dpi": 180, "savefig.bbox": None, "savefig.pad_inches": 0,
    })


def canvas(ncols=2, height=78):
    fig, axes = plt.subplots(1, ncols, figsize=(178/25.4, height/25.4))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.19,
                        top=0.81, wspace=0.36)
    for i, ax in enumerate(axes):
        ax.text(0, 1.08, f"({chr(97+i)})", transform=ax.transAxes,
                fontsize=8, fontweight="bold", va="bottom")
    return fig, axes


def save(fig, outdir, name):
    fig.savefig(outdir / f"{name}.pdf", bbox_inches=None)
    svg = outdir / f"{name}.svg"
    fig.savefig(svg, bbox_inches=None)
    svg.write_text("\n".join(line.rstrip() for line in
                            svg.read_text(encoding="utf-8").splitlines()) + "\n",
                   encoding="utf-8", newline="\n")
    fig.savefig(outdir / f"{name}.png", dpi=600, bbox_inches=None)
    outputs = [str(outdir/f"{name}.{ext}") for ext in ("pdf", "svg", "png")]
    plt.close(fig)
    return outputs


def rows(table, generator, theta=None, L=None):
    return [r for r in table if r["generator"] == generator
            and (theta is None or abs(r["theta"]-theta) < 1e-4)
            and (L is None or r["L"] == L)]


def errorbar(ax, x, y, err, generator):
    label, color, marker = STYLE[generator]
    return ax.errorbar(x, y, yerr=err, color=color, marker=marker, ms=3.1,
                       capsize=1.6, elinewidth=0.75, label=label)


def positive(values, label):
    values = np.asarray(values)
    if not np.all(np.isfinite(values) & (values > 0)):
        raise ValueError(f"{label} must be finite and positive on log axes")


def fig1(noise, outdir):
    """Production covariance, same-sample layout control and exact raw AR sums."""
    fig, axes = canvas(3, 122)
    for i, ax in enumerate(axes):
        ax.set_position([0.08 + i * 0.315, 0.12, 0.245, 0.285])
        ax.texts[0].set_text(f"({chr(99+i)})")
    matched = rows(noise["series"], "fgn", 0.25, 32)[0]
    for i, (key, ylabel) in enumerate((
            ("time", "Temporal correlation C(Δt = 1)"),
            ("space", "Spatial correlation C(Δx = 1)"))):
        ax = fig.add_axes([0.08 + 0.49*i, 0.625, 0.405, 0.30])
        means = [r[f"C_{key}"][0] for r in (matched, noise["legacy"])]
        errors = [r[f"C_{key}_err"][0] for r in (matched, noise["legacy"])]
        bars = ax.bar([0, 1], means, yerr=errors, width=0.52,
                      color=[STYLE["fgn"][1], STYLE["shared"][1]],
                      edgecolor="white", linewidth=0.5, capsize=3,
                      error_kw={"elinewidth": 0.8, "ecolor": NEUTRAL})
        bars[1].set_hatch("///")
        ax.axhline(0, color=NEUTRAL, lw=0.6)
        if key == "time":
            ax.axhline(matched["target_time"][0], color=NEUTRAL, ls="--", lw=0.7)
        for x, mean, se in zip((0, 1), means, errors):
            ax.text(x, mean+se+0.012, f"{mean:.4f} ± {se:.4f}",
                    ha="center", va="bottom", fontsize=7)
        ax.set(ylabel=ylabel, ylim=(-0.015, 0.335), xlim=(-0.6, 1.6),
               xticks=[0, 1], xticklabels=["Correct layout", "Axis-permuted"],
               yticks=[0, 0.1, 0.2, 0.3])
        ax.tick_params(axis="x", length=0, pad=6)
        ax.text(0, 1.04, f"({chr(97+i)})", transform=ax.transAxes,
                fontweight="bold", fontsize=8, va="bottom")
    current = rows(noise["series"], "fgn", 0.25, 64)[0]
    ax = axes[0]
    positive(current["C_time"], "Temporal correlation")
    positive(current["target_time"], "Target temporal correlation")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.errorbar(noise["lags"], current["C_time"], yerr=current["C_time_err"],
                color=STYLE["fgn"][1], marker="^", ms=3, capsize=1.5,
                label="Production")
    ax.plot(noise["lags"], current["target_time"], "--", color=NEUTRAL, label="Target")
    ax.set(xlabel="Temporal lag (sweeps)", ylabel="Binary correlation")
    ax.legend(loc="lower left")
    ax = axes[1]
    matched = rows(noise["series"], "fgn", 0.25, 32)[0]
    legacy = noise["legacy"]
    if legacy["L"] != 32 or abs(legacy["theta"]-0.25) > 1e-4:
        raise ValueError("Layout control requires matched theta and L")
    for row, label, color, marker in (
            (matched, "Current layout", STYLE["fgn"][1], "^"),
            (legacy, "Axis-permuted", STYLE["shared"][1], "o")):
        ax.errorbar(noise["spatial_lags"], row["C_space"], yerr=row["C_space_err"],
                    color=color, marker=marker, ms=3, capsize=1.5, label=label)
    ax.axhline(0, color=NEUTRAL, lw=0.7)
    ax.set(xlabel="Spatial lag (sites)", ylabel="Binary correlation")
    ax.set_xticks(noise["spatial_lags"])
    ax.legend(loc="upper right")
    ax = axes[2]
    ar = noise["ar_zero"]
    for key, label, color, ls in (
            ("C_shared", "Shared innovations", STYLE["shared"][1], "-"),
            ("C_independent", "Independent innovations", NEUTRAL, "--")):
        positive(ar[key], key)
        ax.loglog(ar["lags"], ar[key], color=color, ls=ls, label=label)
    ax.set(xlabel="Temporal lag (sweeps)", ylabel="Raw-field correlation")
    ax.legend(loc="lower left")
    return save(fig, outdir, "fig1_generator")


def fig2(summary, outdir):
    fig, axes = canvas()
    for generator in STYLE:
        data = sorted(rows(summary["exponents"], generator), key=lambda r: r["theta"])
        errorbar(axes[0], [r["theta"] for r in data], [r["alpha"] for r in data],
                 [r["alpha_err"] for r in data], generator)
        pairs = rows(summary["exponents"], generator, 0.25)[0]["alpha_pairs"]
        pairs = sorted(pairs, key=lambda r: r["midpoint"])
        errorbar(axes[1], [r["midpoint"] for r in pairs], [r["alpha"] for r in pairs],
                 [r["alpha_err"] for r in pairs], generator)
    for ax in axes:
        ax.axhline(0.5, color=NEUTRAL, lw=0.8)
    axes[0].set(xlabel="θ", ylabel=r"Finite-size slope $\alpha_{\mathrm{fit}}$")
    axes[1].set(xlabel="Pair midpoint L (sites)", ylabel=r"$\alpha_{\mathrm{eff}}$",
                title="θ = 0.25")
    axes[1].set_xscale("log", base=2)
    axes[1].xaxis.set_major_formatter(ScalarFormatter())
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.53, 1.0))
    return save(fig, outdir, "fig2_exponents")


def fig3(summary, outdir):
    fig, axes = canvas(3, height=80)
    fig.subplots_adjust(left=0.07, right=0.985, wspace=0.50)
    for generator in STYLE:
        history, _, _ = profile_history(summary, generator)
        h = history[:, -1]
        samples = np.abs(np.fft.rfft(h, axis=1))**2/h.shape[1]
        spectrum = samples[:, 1:].mean(axis=0)
        se = samples[:, 1:].std(axis=0, ddof=1)/np.sqrt(h.shape[0])
        positive(spectrum-se, "Spectrum minus one run SE")
        modes = np.arange(1, len(spectrum)+1)
        label, color, _ = STYLE[generator]
        axes[2].loglog(modes, spectrum, color=color, label=label)
        axes[2].fill_between(modes, spectrum-se, spectrum+se, color=color, alpha=0.16)
    axes[2].set(xlabel="Mode m (k = 2πm/L)", ylabel="Ensemble S(k)",
                title="θ = 0.40, L = 256")
    for generator in STYLE:
        data = sorted(rows(summary["cells"], generator, L=256), key=lambda r: r["theta"])
        errorbar(axes[0], [r["theta"] for r in data], [r["p_flat"] for r in data],
                 [r["p_flat_err"] for r in data], generator)
    axes[0].set(xlabel="θ", ylabel=r"$p_{\mathrm{flat}}$", title="L = 256")
    axes[0].legend(loc="best")
    for theta, color, marker in ((0.1, "#8fc1df", "o"),
                                 (0.25, "#3c8fc4", "s"),
                                 (0.4, "#123e7d", "^")):
        data = sorted(rows(summary["cells"], "fgn", theta), key=lambda r: r["L"])
        axes[1].errorbar([r["L"] for r in data], [r["dH"] for r in data],
                         yerr=[r["dH_err"] for r in data], color=color,
                         marker=marker, ms=3.1, capsize=1.6,
                         elinewidth=0.75, label=f"θ = {theta:g}")
    axes[1].set_xscale("log", base=2)
    axes[1].xaxis.set_major_formatter(ScalarFormatter())
    axes[1].axhline(0, color=NEUTRAL, lw=0.7)
    axes[1].set(xlabel="L (lattice sites)", ylabel="Local-pmf excess ΔH", title="fGn spectrum")
    axes[1].legend(loc="best")
    return save(fig, outdir, "fig3_fixed_theta")


def profile_history(summary, generator, theta=0.4):
    cell = rows(summary["cells"], generator, theta, 256)[0]
    path = (ROOT / cell["raw_file"]).resolve()
    expected = (ROOT / "data/raw" / summary["protocol"]["tags"][generator]).resolve()
    if path.parent != expected:
        raise ValueError(f"Raw file is outside the declared {generator} campaign")
    with np.load(path, allow_pickle=False) as raw:
        if generator != "shared" and str(raw["noise_layout"]) != "spectral-run-time-site-v2":
            raise ValueError("Morphology requires corrected-layout spectral data")
        h = raw["profiles"].astype(float)
        times = raw["check_times"].copy()
        spec = json.loads(str(raw["spec"]))
    if h.shape != (summary["protocol"]["n_runs"], 32, 256):
        raise ValueError("Morphology must include all declared independent runs")
    assert np.all(np.diff(times) > 0) and times[0] > 0
    h -= h.mean(axis=2, keepdims=True)
    widths = np.sqrt(np.mean(h[:, -1]**2, axis=1))
    index = int(np.argmin(np.abs(widths-np.median(widths))))
    source = {"generator": generator, "theta": theta, "L": 256,
              "raw_file": cell["raw_file"],
              "raw_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "representative_run_index": index,
              "selection": "final width nearest the all-run median; first tie",
              "quantitative_n_runs": int(h.shape[0]),
              "displayed_snapshots": len(times), "displayed_sites": h.shape[2],
              "T": spec["T"], "check_times": times.tolist()}
    return h, times, source


def surface_panel(fig, rect, history, times, source, norm, panel):
    """Piecewise continuous surface through every stored x/time height sample."""
    ax = fig.add_axes(rect, projection="3d")
    x, logt = np.meshgrid(np.arange(256), np.log10(times / source["T"]))
    height = history[source["representative_run_index"]]
    colors = LightSource(azdeg=315, altdeg=42).shade(
        height, cmap=mpl.colormaps["viridis"], norm=norm,
        vert_exag=0.45, blend_mode="soft")
    ax.plot_surface(x, logt, height, rcount=len(times), ccount=256,
                    facecolors=colors, linewidth=0, antialiased=False,
                    shade=False, rasterized=True)
    ax.set(xlim=(0, 255), ylim=(logt.min(), 0), zlim=(norm.vmin, norm.vmax),
           xticks=[0, 128, 255], yticks=[-3, -1, 0], zticks=[-60, 0, 60],
           xlabel=r"$x$", ylabel=r"$\log_{10}(t/T)$", zlabel=r"$h-\bar h$")
    ax.view_init(elev=32, azim=-62)
    ax.set_box_aspect((1.65, 1.35, 0.95), zoom=0.92)
    ax.tick_params(labelsize=6.5, pad=0, width=0.5)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.label.set_size(8)
        axis.labelpad = 1
        axis.line.set_linewidth(0.5)
        axis.set_pane_color((0.97, 0.98, 0.985, 0.7))
        axis._axinfo["grid"].update(color=(0.77, 0.82, 0.85, 0.6), linewidth=0.35)
    ax.text2D(0.0, 0.98, f"({panel})", transform=ax.transAxes,
              fontsize=8, fontweight="bold", va="top")


def fig4(summary, outdir):
    """Three continuous surfaces with coordinates and panel labels, no titles."""
    histories = {g: profile_history(summary, g) for g in STYLE}
    bound = max(np.abs(h[s["representative_run_index"]]).max()
                for h, _, s in histories.values())
    bound = 10 * np.ceil(bound/10)
    norm = Normalize(-bound, bound)
    fig = plt.figure(figsize=(178/25.4, 75/25.4))
    selected = []
    for column, generator in enumerate(STYLE):
        history, times, source = histories[generator]
        surface_panel(fig, [0.005+column*0.326, 0.17, 0.275, 0.73],
                      history, times, source, norm, chr(97+column))
        selected.append(source)
    return save(fig, outdir, "fig4_morphology"), selected


def fig6(summary, outdir):
    """Nine labelled surfaces, ordered by theta (rows) and input (columns)."""
    histories = {(g, theta): profile_history(summary, g, theta)
                 for theta in (0.1, 0.25, 0.4) for g in STYLE}
    bound = max(np.abs(h[s["representative_run_index"]]).max()
                for h, _, s in histories.values())
    bound = 10 * np.ceil(bound/10)
    norm = Normalize(-bound, bound)
    fig = plt.figure(figsize=(178/25.4, 165/25.4))
    selected = []
    for row, theta in enumerate((0.1, 0.25, 0.4)):
        for col, generator in enumerate(STYLE):
            history, times, source = histories[generator, theta]
            surface_panel(fig, [0.025+col*0.326, 0.715-row*0.318, 0.275, 0.25],
                          history, times, source, norm, chr(97+row*3+col))
            selected.append(source)
    return save(fig, outdir, "fig6_morphology_atlas"), selected


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path,
                        default=ROOT/"data/analysis/revision_20261026/summary.json")
    parser.add_argument("--outdir", type=Path, default=ROOT/"figs/revision_20261026")
    args = parser.parse_args(argv)
    summary_path = args.summary.resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    configure_style()
    noise_path = summary_path.with_name("noise_validation.json")
    noise = json.loads(noise_path.read_text(encoding="utf-8"))
    outputs = fig1(noise, outdir) + fig2(summary, outdir) + fig3(summary, outdir)
    morphology, selected = fig4(summary, outdir)
    outputs += morphology
    atlas, atlas_selected = fig6(summary, outdir)
    outputs += atlas
    sources = {"summary": summary_path.relative_to(ROOT).as_posix() if summary_path.is_relative_to(ROOT) else str(summary_path),
               "summary_sha256": hashlib.sha256(summary_path.read_text(encoding="utf-8").encode("utf-8")).hexdigest(),
               "noise_validation": noise_path.relative_to(ROOT).as_posix() if noise_path.is_relative_to(ROOT) else str(noise_path),
               "noise_sha256": hashlib.sha256(noise_path.read_text(encoding="utf-8").encode("utf-8")).hexdigest(),
               "protocol": summary["protocol"], "morphology_sources": selected,
               "atlas_sources": atlas_selected,
               "morphology_rendering": {
                   "type": "continuous shaded surfaces",
                   "coordinates": ["x", "log10(t/T)", "h - mean(h)"],
                   "height_limits": [-60, 60], "colormap": "viridis",
                   "vertices": "all 32 stored times and all 256 sites; linear faces",
                   "annotations": "coordinates, ticks and panel letters; no titles",
                   "quantitative_spectrum": "panel c of fig3_fixed_theta"},
               "noise_comparison": {"series": "fgn theta=0.25 L=32",
                                    "n_runs": 8, "lag": 1,
                                    "panels": "temporal and spatial mean +/- one run SE",
                                    "negative_control": "axis-permuted same samples"},
               "outputs": [Path(p).relative_to(ROOT).as_posix() if Path(p).is_relative_to(ROOT) else p for p in outputs],
               "uncertainty": "Supplied run errors; S(k) uses all-run mean and SE.",
               "text_hash_normalization": "UTF-8 with LF line endings; binary archive hashes remain byte-exact."}
    (outdir/"figure_sources.json").write_text(json.dumps(sources, indent=2), encoding="utf-8")
    print(f"Exported {len(outputs)} files at 178 mm width to {outdir}")


if __name__ == "__main__":
    main()
