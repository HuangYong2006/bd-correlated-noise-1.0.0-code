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
    fig, axes = canvas(3, 80)
    fig.subplots_adjust(left=0.08, wspace=0.43)
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
    ax.set(xlabel="Temporal lag (sweeps)", ylabel="Binary correlation",
           title="fGn spectrum, θ = 0.25, L = 64")
    ax.legend(loc="lower left")
    ax = axes[1]
    matched = rows(noise["series"], "fgn", 0.25, 32)[0]
    legacy = noise["legacy"]
    if legacy["L"] != 32 or abs(legacy["theta"]-0.25) > 1e-4:
        raise ValueError("Layout control requires matched theta and L")
    for row, label, color, marker in (
            (matched, "Current layout", STYLE["fgn"][1], "^"),
            (legacy, "Legacy layout", STYLE["shared"][1], "o")):
        ax.errorbar(noise["spatial_lags"], row["C_space"], yerr=row["C_space_err"],
                    color=color, marker=marker, ms=3, capsize=1.5, label=label)
    ax.axhline(0, color=NEUTRAL, lw=0.7)
    ax.set(xlabel="Spatial lag (sites)", ylabel="Binary correlation",
           title="Same samples, θ = 0.25, L = 32")
    ax.set_xticks(noise["spatial_lags"])
    ax.legend(loc="upper right")
    ax = axes[2]
    ar = noise["ar_zero"]
    for key, label, color, ls in (
            ("C_shared", "Shared innovations", STYLE["shared"][1], "-"),
            ("C_independent", "Independent innovations", NEUTRAL, "--")):
        positive(ar[key], key)
        ax.loglog(ar["lags"], ar[key], color=color, ls=ls, label=label)
    ax.set(xlabel="Temporal lag (sweeps)", ylabel="Raw-field correlation",
           title="Exact AR covariance, θ = 0")
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
    fig, axes = canvas()
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


def profiles(summary, generator):
    cell = rows(summary["cells"], generator, 0.4, 256)[0]
    path = (ROOT / cell["raw_file"]).resolve()
    expected = (ROOT / "data/raw" / summary["protocol"]["tags"][generator]).resolve()
    if path.parent != expected:
        raise ValueError(f"Raw file is outside the declared {generator} campaign")
    with np.load(path, allow_pickle=False) as raw:
        if generator != "shared" and str(raw["noise_layout"]) != "spectral-run-time-site-v2":
            raise ValueError("Morphology requires corrected-layout spectral data")
        h = raw["profiles"][:, -1, :].astype(float)
    if h.shape != (summary["protocol"]["n_runs"], 256):
        raise ValueError("Morphology must include all declared independent runs")
    h -= h.mean(axis=1, keepdims=True)
    return h, cell["raw_file"]


def fig4(summary, outdir):
    fig, axes = canvas(height=84)
    offset = None
    selected = []
    for generator in STYLE:
        label, color, _ = STYLE[generator]
        h_all, raw_file = profiles(summary, generator)
        widths = np.sqrt(np.mean(h_all**2, axis=1))
        index = int(np.argmin(np.abs(widths-np.median(widths))))
        h = h_all[index]
        offset = 0.0 if offset is None else bottom - 7.0 - float(h.max())
        axes[0].plot(np.arange(h.size), h+offset, color=color, lw=0.7, label=label)
        axes[0].axhline(offset, color="#b3b3b3", ls=":", lw=0.5)
        bottom = float(h.min()+offset)
        selected.append({"generator": generator, "raw_file": raw_file,
                         "representative_run_index": index,
                         "selection": "final width nearest the all-run median",
                         "quantitative_n_runs": int(h_all.shape[0])})
        spectra = np.abs(np.fft.rfft(h_all, axis=1))**2/h_all.shape[1]
        spectrum = spectra[:, 1:].mean(axis=0)
        error = spectra[:, 1:].std(axis=0, ddof=1)/np.sqrt(h_all.shape[0])
        positive(spectrum-error, "Spectrum minus one run SE")
        m = np.arange(1, spectrum.size+1)
        axes[1].loglog(m, spectrum, color=color, lw=1.0, label=label)
        axes[1].fill_between(m, spectrum-error, spectrum+error, color=color, alpha=0.16)
    axes[0].set(xlabel="x (lattice sites)", ylabel="Mean-centred height (offset)",
                title="θ = 0.40, L = 256", xlim=(0, 255))
    axes[0].set_yticks([])
    top = max(np.max(line.get_ydata()) for line in axes[0].lines)
    axes[0].plot([242, 242], [top+5, top+15], color=NEUTRAL, lw=0.8)
    axes[0].text(236, top+10, "10", ha="right", va="center", fontsize=6.2)
    axes[0].set_ylim(bottom-5, top+20)
    axes[1].set(xlabel="Mode m (k = 2πm/L)", ylabel="Ensemble S(k)",
                title="All-run mean ± one SE")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.53, 1.0))
    return save(fig, outdir, "fig4_morphology"), selected


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
    sources = {"summary": summary_path.relative_to(ROOT).as_posix() if summary_path.is_relative_to(ROOT) else str(summary_path),
               "summary_sha256": hashlib.sha256(summary_path.read_text(encoding="utf-8").encode("utf-8")).hexdigest(),
               "noise_validation": noise_path.relative_to(ROOT).as_posix() if noise_path.is_relative_to(ROOT) else str(noise_path),
               "noise_sha256": hashlib.sha256(noise_path.read_text(encoding="utf-8").encode("utf-8")).hexdigest(),
               "protocol": summary["protocol"], "morphology_sources": selected,
               "outputs": [Path(p).relative_to(ROOT).as_posix() if Path(p).is_relative_to(ROOT) else p for p in outputs],
               "uncertainty": "Supplied run errors; S(k) uses all-run mean and SE.",
               "text_hash_normalization": "UTF-8 with LF line endings; binary archive hashes remain byte-exact."}
    (outdir/"figure_sources.json").write_text(json.dumps(sources, indent=2), encoding="utf-8")
    print(f"Exported {len(outputs)} files at 178 mm width to {outdir}")


if __name__ == "__main__":
    main()
