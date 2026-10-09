# Noise-generator bias in temporally correlated ballistic deposition

[![DOI](https://img.shields.io/badge/DOI-to%20be%20assigned-lightgrey.svg)](#citation)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Data: NPZ](https://img.shields.io/badge/data-NPZ%20archives-informational.svg)](docs/DATA.md)

Pure-numerical study of the (1+1)-dimensional **ballistic deposition (BD)** model driven by
noise that is power-law correlated in time. The work asks whether the critical temporal
correlation exponent `theta_c`, where Kardar–Parisi–Zhang (KPZ) scaling gives way to a
long-range regime, can be located at accessible system sizes — and reports that it cannot,
because the noise generator used by the published literature on this problem does not
sample the short-range regime at all.

> **Scope of the claim.** This repository does **not** report a measured `theta_c`. Within
> `L = 32 … 512`, no crossover is resolved for `theta <= 0.25` with a generator of correct
> correlation amplitude, and the renormalisation-group prediction `theta_th = 1/6` is
> neither confirmed nor contradicted. What is established is that the generator bias
> accounts for the reported thresholds, and that the standard fixed-significance onset
> construction cannot locate a crossover when the detection floor depends on system size.

---

## Contents

- [Key results](#key-results)
- [Model and conventions](#model-and-conventions)
- [Repository layout](#repository-layout)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Reproducing the paper](#reproducing-the-paper)
- [Data](#data)
- [Figures](#figures)
- [Conventions and gotchas](#conventions-and-gotchas)
- [Citation](#citation)
- [License](#license)
- [References](#references)

---

## Key results

### 1. The FFGN generator is discontinuous at `theta = 0`

Every published study on this problem generates the noise with Mandelbrot's **fast
fractional Gaussian noise (FFGN)**. Its normalised correlation amplitude is essentially
`theta`-independent, `A_FFGN ≈ 1.27`, whereas the exact fractional-Gaussian amplitude
`A_fGn = 2 theta (theta + 1/2)` **vanishes** at `theta = 0`. So `FFGN(theta -> 0+)` is not
white noise: we measure `C(1)/C(0) = 0.195` and a `C(tau) ∝ tau^-1` tail, where white noise
requires `C(tau) = 0` for `tau >= 1`.

The effect is **intrinsic, not a truncation artefact** — the fitted amplitude is
`1.3558 / 1.3559 / 1.3559` for `N = 23 / 32 / 48` AR components.

Because the model binarises the noise, the calibration has to be done *after* thresholding
(`C_bin = (2/pi) arcsin C`). Matching the binarised correlation against the binarised exact
fGn family gives:

| `theta_nom` | FFGN `theta_eff` | spectral `theta_eff` |
|---|---|---|
| **0.000** | **0.1888** | **0.0000** |
| 0.050 | 0.2037 | 0.0435 |
| 0.100 | 0.2215 | 0.0887 |
| **1/6** | **0.2492** | **0.1525** |
| 0.250 | 0.2934 | 0.2357 |
| 0.400 | 0.3966 | 0.3825 |

**For every `theta_nom >= 0` the FFGN family lies above the RG threshold `theta_th = 1/6`.**
It never samples the short-range KPZ regime.

Consequences measured *inside the model* (`L = 256`):

| generator | `theta` | `p_flat` | `alpha_s` | `W_sat` |
|---|---|---|---|---|
| i.i.d. (true white) | 0 | 0.4330 | 0.514 | 3.65 |
| FFGN | `theta -> 0` | 0.3359 | — | 5.61 |
| spectral | 0 | 0.4339 | 0.516 | 3.62 |
| spectral | 1/6 | 0.4505 | 0.536 | 3.81 |
| spectral | 0.40 | 0.5110 | 0.548 | 4.05 |
| FFGN | 0.40 | **0.1106** | **1.049** | **17.99** |

FFGN produces a **23 % discontinuity** in the structural order parameter at `theta = 0`,
*inverts* its `theta`-dependence, drives `alpha_s` to 1.05 (reproducing the published
faceting trend) and inflates `W_sat` by a factor of 4.4.

### 2. The model is single-step (structural result)

With 0/1 noise, the BD rule preserves `|Delta h| <= 1` from a flat initial condition (proved
in `manuscript/main.tex`, App. B, and verified in every cell by `code/check_single_step.py`).
Hence the increment density has support on `{-1, 0, +1}`. Moreover `sum_x Delta h(x) = 0`
telescopes exactly on a periodic lattice, so `p_up = p_down` **per configuration** and
`skew(Delta h) ≡ 0` to machine precision — the pmf has exactly one free parameter,
`p_flat`. "Faceting" in this model therefore cannot be overhang-based.

### 3. The standard `theta_c` construction is confounded (a self-correction)

In nominal FFGN units both diagnostics give `theta_c = 0.0002 ± 0.0007` and
`0.0001 ± 0.0011`, but that onset *is* the generator discontinuity (the `theta = 0`
reference is i.i.d. white while every `theta > 0` cell uses FFGN), so it is not a
measurement of `theta_c`.

On the corrected axis, with run counts matched at 96 everywhere, the fixed-significance
onset decreases monotonically over a factor of 16 in `L`:

| `L` | 32 | 64 | 128 | 256 | 512 |
|---|---|---|---|---|---|
| `theta*(z>=3)` | 0.2415 | 0.2336 | 0.1710 | 0.1375 | 0.0823 |
| `theta*(z>=5)` | 0.3306 | 0.3073 | 0.2686 | 0.2410 | 0.1471 |

and a forced-zero finite-size-scaling fit is formally excellent (`chi2/dof = 0.75/3`,
`p = 0.86`). **This drift is an artefact.** The split-half null floor falls from 0.0109 at
`L = 32` to 0.0024 at `L = 512` (more columns per run), so a fixed-`z` crossing moves to
smaller `theta` even for a scale-free observable. Matching run counts is necessary but not
sufficient: the `L`-dependence enters through the *number of columns*.

### 4. The decisive test is at fixed `theta`

A genuine crossover requires the deviation to *grow* with `L`, since larger `L` resolves
more of the long-range behaviour. Over `L = 32 … 512`, `DeltaH(theta, L)` at fixed `theta`
is **`L`-independent within errors** for `theta <= 0.25` (`dDeltaH/dlnL` consistent with
zero for 22 of 26 values of `theta`, no systematic sign), but becomes **significantly
positive for `theta >= 0.30`**. The crossover is real and simply lies beyond `L ≈ 512` for
the smaller values.

### 5. Exponents, with honest error bars

Residual-based errors on a 4–5 point log-log fit understated the uncertainty by ~8×
(`±0.0016` vs `±0.0127` from a run bootstrap). With bootstrap errors:

| `theta` | −0.2 | 0 | 1/6 | 0.25 | 0.40 |
|---|---|---|---|---|---|
| `alpha` (spectral) | 0.4947(14) | 0.5062(13) | 0.5082(12) | 0.5311(12) | 0.5750(12) |
| `alpha` (FFGN) | — | 0.4986(13) | 0.6330(15) | 0.7019(13) | 0.8585(12) |

A supplementary campaign with the **exact fGn spectral density** closed the one systematic
that could have faked a flat `alpha` — the bare power-law filter's amplitude deficit at
large `theta`:

| `theta` | FFGN (published) | power-law filter | exact fGn filter |
|---|---|---|---|
| 0 | 0.4986(13) | 0.5062(13) | 0.5076(12) |
| 1/6 | 0.6330(15) | 0.5082(12) | 0.5089(11) |
| **0.20** | 0.6855(14) | 0.5269(12) | **0.5368(11)** |
| **0.25** | 0.7019(13) | 0.5311(12) | **0.5383(12)** |
| 0.40 | 0.8585(12) | 0.5750(12) | 0.5724(12) |

The two corrected filters agree to within **0.008** everywhere, while FFGN lies **0.16–0.29**
higher. `alpha_eff` from consecutive size pairs decreases towards 1/2 as `L` grows for the
spectral generator at `theta >= 1/6`, but not under FFGN — the corrected deviations behave
like finite-size corrections; the FFGN ones behave asymptotically.

**Negative-`theta` control** (never simulated before): `alpha = 0.495 ± 0.014`,
`0.505 ± 0.012`, `0.499 ± 0.009` at `theta = −0.2, −0.1, −0.05`, all consistent with 1/2 as
theory requires, while `p_flat` does drop — separating "the noise changed" from "the
universality class changed".

### 6. Diagnostics: three of four candidate observables were rejected *with data*

**Diagnostic 1 — Hellinger distance (distributional).** `H² = 1 − Σ sqrt(p_k q_k)` between
the **exact discrete pmf** of the nearest-neighbour increment, with a split-half null floor
measured at matched ensemble size.

| observable | null floor | verdict |
|---|---|---|
| **increment pmf (exact, binning-free)** | **0.0024–0.0109** | accepted, `z` up to 30 |
| terrace-length pmf (exact) | 0.024–0.099 | usable, noisier |
| histogrammed standardised `Delta h_r` | 0.080–0.181 | rejected, signal ≈ 0 |
| histogrammed standardised height | 0.102–0.251 | rejected, blind control |

Two traps worth recording: **histogramming a discrete distribution before computing a
Hellinger distance is invalid** (spikes land in disjoint bins, `H` saturates near its
maximum regardless of the physical difference), and **any observable retaining the
fluctuation amplitude reports a transition at `theta = 0+`**.

**Diagnostic 2 — Isolation Forest (multivariate geometry).** Trained on the `theta = 0`
ensemble of amplitude-free scaling/morphology features, threshold calibrated at a fixed
false-positive rate on held-out `theta = 0` data. Its features are exactly the
roughness-scaling quantities from which `alpha` is built, so in the corrected track it
correctly reports *no* restructuring (`|z| < 2` for all `theta <= 0.4`), while the
distributional diagnostic responds strongly. That contrast is itself a result: over the
accessible range the change in the interface is a local one-parameter effect, not a
multi-scale one.

---

## Model and conventions

Noise correlator (Squizzato & Canet 2019, Eq. 4), equivalently a power-law temporal
correlation:

```
D(omega, k) = D0 + D_theta * omega^(-2 theta)   <=>   <eta(x,t) eta(x',t')> = 2 D delta_{x,x'} |t - t'|^(2 theta - 1)
```

Discrete model (Alés & López 2019; Wang & Xia 2025; Jiang & Xia 2026), parallel update, flat
initial condition, periodic boundaries:

```
h(x, t+1) = max[ h(x,t) + eta'_theta(x,t),  h(x-1,t),  h(x+1,t) ]
eta'_theta(x,t) = 1 if eta_theta(x,t) > 0 else 0,     sites independent
```

| symbol | meaning |
|---|---|
| `theta` | nominal temporal-correlation exponent of the noise |
| `theta_eff` | effective exponent after binarisation, matched to exact fGn |
| `theta_th = 1/6` | RG threshold between KPZ and long-range scaling |
| `L` | number of columns (system size) |
| `T` | number of sweeps; `T = t_mult * L^1.5` |
| `W_sat` | saturated interface width |
| `alpha`, `beta` | roughness and growth exponents |
| `p_flat` | structural order parameter: fraction of `Delta h = 0` |

Two noise generators are provided, because the effective exponent at a given nominal
`theta` depends on the generator:

1. **`ffgn`** — Mandelbrot's FFGN, the generator used by the published literature.
   Streaming, `O(N)` state per site, arbitrarily long runs are cheap in memory. This is the
   default and the literature-comparable choice. **Discontinuous at `theta = 0`.**
2. **`spectral`** — exact spectral filtering plus Gaussian thresholding, continuous at
   `theta = 0`. Two filter amplitudes are available:
   * `powerlaw` — `|f|^(-theta)`; the family used for the main campaign. `theta` axis
     accurate to ~10 %, under-correlates at large `theta`.
   * `fgn` — the exact fractional-Gaussian spectral density; `theta` axis accurate to ~1 %,
     exactly white at `theta = 0`.

---

## Repository layout

```
.
├── code/                        simulation, diagnostics, analysis (flat import namespace)
│   ├── bd_core.py               model, FFGN + 2 spectral generators, numba simulator
│   ├── run_campaign.py          campaign runner: writes one .npz per (theta, L) cell
│   ├── features.py              amplitude-free dimensionless features
│   ├── diagnostics.py           Hellinger (exact pmf, null-calibrated) + Isolation Forest
│   ├── analyze.py               driver: diagnostics, onsets, FSS fits, exponents
│   ├── fixed_theta_scaling.py   the decisive test (fixed-theta L-independence)
│   ├── alpha_drift.py           alpha_eff(L) from size pairs + bootstrap errors
│   ├── feature_ablation.py      diagnostic correlation + feature-family ablation
│   ├── calibrate_binarised.py   the generator calibration (binarised, both families)
│   ├── calibrate_spectral.py    calibration of the power-law spectral family
│   ├── calibrate_generator.py   first-generation Gaussian-domain calibration (superseded)
│   ├── validate_fgn_spectrum.py exact-fGn-spectrum filter validation
│   ├── validate_model.py        model + FFGN validation  (run this first)
│   ├── compare_observables.py   head-to-head test of candidate observables
│   ├── check_single_step.py     verification of |Delta h| <= 1 and skew ≡ 0
│   ├── check_generator_amplitude.py  A_FFGN vs A_fGn, N-convergence
│   ├── check_generator_limit.py FFGN converged in the number of components
│   ├── quick_hellinger_scan.py  fast scan of the distributional diagnostic
│   ├── bench_spectral.py        timing benchmark for the spectral generator
│   ├── figures_ctp.py           the manuscript figures (4, at 178 mm)
│   └── figures.py, fig_*.py     legacy figure scripts, superseded by figures_ctp.py
├── data/
│   ├── raw/{main,spec,fgnspec}/ one .npz per (theta, L) + manifest_bd.json
│   ├── raw/{explore,testspec}/  preliminary scouting campaigns
│   ├── analysis/<tag>_bd/       diagnostics.csv, summary.json, rows.npz
│   └── *.txt, *.npz             calibration and validation outputs
├── figs/
│   ├── ctp/                     the four manuscript figures (pdf/svg/tiff + audit json)
│   └── {generator,physics}/     legacy figures, superseded, not included in the paper
├── manuscript/
│   ├── main.tex, refs.bib       the paper (REVTeX 4.2, PRE two-column)
│   ├── abstract_zh.md           Chinese abstract and notes
│   ├── cover_letter.md          submission cover letter
│   ├── overleaf/                Overleaf-ready build tree
│   └── overleaf_bd_paper.zip    upload this to Overleaf
├── docs/                        REPRODUCE.md, DATA.md, CODE_MAP.md, LEGACY.md, review notes
├── tests/                       unit tests for the core model invariants
├── .github/workflows/tests.yml  CI: byte-compile + unit tests
├── pack_repo.ps1                build the upload archive
├── CITATION.cff, LICENSE, CHANGELOG.md
├── requirements.txt, pyproject.toml
├── README.md                    this file (English)
└── README_zh.md                 Chinese companion overview
```

A per-file catalogue, including which scripts are live and which are superseded, is in
[`docs/CODE_MAP.md`](docs/CODE_MAP.md), with the rationale for the retained dead scripts in
[`docs/LEGACY.md`](docs/LEGACY.md). The statistical pre-submission review notes kept at the
repository root are indexed in [`docs/REVIEW_NOTES.md`](docs/REVIEW_NOTES.md).

---

## Installation

Python **3.10+** is required. The pipeline is CPU-only; no GPU is used.

```bash
git clone https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code.git
cd bd-correlated-noise-1.0.0-code
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux / macOS
source .venv/bin/activate

pip install -r requirements.txt
```

Core dependencies: `numpy`, `scipy`, `numba` (JIT + `prange` threads), `scikit-learn`
(Isolation Forest), `matplotlib` (figures).

Optional:

- `PyMuPDF` (`fitz`) — only for the PDF-inspection utilities in `tools/`.
- `pytest` — only for the test suite in `tests/`.

Verify the install and the model:

```bash
python code/validate_model.py
```

This checks the noise correlator, reproduces `alpha = 1/2` and `beta = 1/3` at `theta = 0`,
and reports timings used to plan the campaign cost.

---

## Quick start

Run a small campaign and analyse it:

```bash
# 1. simulate: 1 size x 3 exponents x 64 runs, a few minutes on a laptop
python code/run_campaign.py --tag demo --rule bd --L 64 \
    --thetas 0 0.1667 0.40 --runs 64 --t-mult 15

# 2. analyse: onsets, finite-size scaling fits, exponents
python code/analyze.py --tag demo --rule bd --checkpoints 4

# 3. inspect the result tables
#    data/analysis/demo_bd/{diagnostics.csv, summary.json, rows.npz}
```

Note that onset extraction and finite-size scaling need **at least three sizes**, so the
`demo` run above only exercises the code path; use the sizes in
[Reproducing the paper](#reproducing-the-paper) for physics.

Model parameters are `SimSpec` fields in `code/bd_core.py`: `theta`, `L`, `T`, `n_runs`,
`rule`, `n_check`, `n_w`, `seed`, `spectrum`.

---

## Reproducing the paper

Full step-by-step instructions, including expected runtimes and the exact configuration of
every campaign, are in [`docs/REPRODUCE.md`](docs/REPRODUCE.md). The short version:

```bash
# --- validation (minutes) -------------------------------------------------
python code/check_single_step.py            # |Delta h| <= 1 and skew ≡ 0
python code/check_generator_amplitude.py    # the FFGN amplitude defect
python code/calibrate_binarised.py          # theta_eff of both generators
python code/validate_fgn_spectrum.py        # exact-fGn-spectrum filter
python code/validate_model.py               # alpha = 1/2, beta = 1/3 at theta = 0
python code/compare_observables.py --tag spec --L 128   # which observable works

# --- campaigns (~5 h on 28 cores total) -----------------------------------
python code/run_campaign.py --tag main  --rule bd --L 64 128 256 --runs 256 \
    --t-mult 15 --thetas 0 0.0025 ... 0.40
python code/run_campaign.py --tag spec  --generator spectral --rule bd \
    --L 32 64 128 256 --runs 192 --t-mult 15 --thetas ...

# --- analysis and figures -------------------------------------------------
python code/analyze.py --tag spec --rule bd --checkpoints 4 --if-runs 96 --hell-runs 96
python code/fixed_theta_scaling.py
python code/alpha_drift.py
python code/figures_ctp.py

# --- manuscript pre-flight, then the upload package -----------------------
python tools/check_tex.py manuscript/main.tex
python tools/table_width.py manuscript/main.tex
python tools/build_overleaf.py
```

Numerical results are deterministic given `--seed`: the simulator uses a counter-based
splitmix64 RNG indexed by `(run, t, x)`, so output does not depend on the number of
threads.

---

## Data

Simulation output is stored as compressed `.npz` archives, one file per `(theta, L)` cell,
named `bd_p0d1667_L512.npz` (sign encoded as `p`/`m`, decimal point as `d`). Each archive
contains the height profiles at a geometric ladder of checkpoints, the interface-width time
series, the per-run statistics and the full `SimSpec` used to produce it.

| tag | generator | `theta` | `L` | runs | purpose |
|---|---|---|---|---|---|
| `main` | FFGN | 0 … 0.40 (22 values) | 64 … 512 | 256 | literature-comparable axis |
| `spec` | spectral, power-law | 0 … 0.40 (26 values) | 32 … 512 | 192 | corrected axis |
| `fgnspec` | spectral, exact fGn | 0 … 0.40 (8 values) | 32 … 512 | 40 cells | systematic closure |
| `explore`, `testspec` | mixed | scouting | 64 … 1024 | — | preliminary scans |

The format, the field-by-field schema, the checksums and the provenance of every cell are
documented in [`docs/DATA.md`](docs/DATA.md). Raw data are **not** tracked by git; see that
document for how the archives are versioned and regenerated.

---

## Figures

`code/figures_ctp.py` renders all four manuscript figures at **exactly 178 mm wide**, which
is within 1 % of the full REVTeX two-column text width (508 pt), so they are placed at
natural size as full-width `figure*` floats with a 7 pt base font and editable TrueType text
(`pdf.fonttype=42`). An invisible full-canvas rectangle pins the `bbox_inches` tight box so
every page comes out at the same width — without it the four figures differed by up to 7 mm
and their fonts would have been scaled differently by `\includegraphics[width=\textwidth]`.

```bash
python code/figures_ctp.py                     # writes figs/ctp/*.pdf/.svg/.tiff
```

Every figure passed a rendered-geometry **collision audit** (0 blocking findings) before the
manuscript package was rebuilt. The audit catches the failures that are invisible in a diff
and easy to miss on screen: text crossed by a data stroke, legends crossed by curves, labels
colliding with tick labels, and content clipped at the canvas edge. The audit helpers live
in an external tool directory; if it is absent the gate reports itself as skipped and the
export still proceeds. Point `DSH_FIGURE_SKILL` at your own copy to re-enable it.

---

## Conventions and gotchas

Hard-won constraints that are easy to reintroduce. Each one cost a debugging session.

- **Parallelism uses numba `prange` threads and threaded `scipy.fft` only.** Process-based
  parallelism (`joblib`/`loky`) interacts badly with confined sandboxes and with the
  counter-based RNG, so do not reintroduce it.
- **Use 5-smooth FFT lengths.** `T = 21722 = 2 x prime` is ~10x slower in `pocketfft`;
  `bd_core.smooth_ceil` returns the smallest 5-smooth integer `>= n`.
- **Cost scaling.** At `theta = 0` the model gives `alpha = 0.4964` and `beta -> 1/3`,
  saturating at `t_x ≈ 0.6 L^1.5` sweeps, so campaign cost scales as `L^2.5` per run.
- **The RSOS rule under a parallel update freezes** and is unusable as-is; BD is the model
  the prior literature uses. `RULE_RSOS` is retained only for completeness.
- **Never histogram a discrete distribution before a Hellinger distance** — the spikes land
  in disjoint bins and `H` saturates. Use the exact pmf.
- **Never match run counts without matching column counts** when comparing a
  significance-based onset across system sizes: the null floor scales with the number of
  columns per run, which is `L` itself.
- **`pdflatex` is not required** to validate the manuscript: `tools/check_tex.py`
  (environments, cross-refs, citations, uncited floats, figure resolution, math/brace
  balance) and `tools/table_width.py` (tabular widths against the 246 pt PRE column) are
  the pre-flight gates used here, plus the rendered collision audit on every figure.
- **Do not edit the `.tex` or `.npz` payloads with `Get-Content`/`Set-Content` round-trips
  in PowerShell** — they mangle UTF-8 and can corrupt compressed archives. Use proper file
  tools.

---

## Citation

If you use this code or data, please cite the manuscript and the archived release. Machine-
readable metadata is in [`CITATION.cff`](CITATION.cff).

```bibtex
@article{Huang2026noisebias,
  title   = {Noise-generator bias in temporally correlated ballistic-deposition simulations},
  author  = {Huang, Yong},
  journal = {Physical Review E},
  year    = {2026},
  note    = {submitted}
}
```

BibTeX entries for the work this study builds on are collected in
`manuscript/refs.bib`.

---

## License

Released under the [MIT License](LICENSE). Copyright (c) 2026 Yong Huang.

---

## References

Key literature this study engages with:

- Medina, Hwa, Kardar & Zhang, *Phys. Rev. A* **39**, 3053 (1989)
- Lam, Sander & Wolf, *Phys. Rev. A* **46**, R6128 (1992)
- Ma & Ma, *Phys. Rev. E* **47**, 3738 (1993)
- Katzav & Schwartz, *Phys. Rev. E* **70**, 011601 (2004)
- Fedorenko, *Phys. Rev. B* **77**, 094203 (2008)
- Alés & López, *Phys. Rev. E* **99**, 062139 (2019)
- Squizzato & Canet, *Phys. Rev. E* **100**, 062143 (2019)
- Song & Xia, *Phys. Rev. E* **103**, 012121 (2021)
- Song & Xia, *J. Stat. Mech.* (2021) 073203
- Hu, Hao & Xia, *Physica A* **619**, 128744 (2023)
- Rodríguez-Fernández et al., *Phys. Rev. E* **110**, 024104 (2024)
- Wang & Xia, *Phys. Rev. E* **111**, 024124 (2025)
- Jiang & Xia, *Commun. Theor. Phys.* **78**, 105602 (2026)
