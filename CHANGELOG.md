# Changelog

All notable changes to this repository are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Repository packaging for public release: `README.md` rewritten as a full research-repository
  entry point; `LICENSE` (MIT), `CITATION.cff`, `requirements.txt`, `pyproject.toml`,
  `.gitignore`, `.gitattributes`, `CHANGELOG.md`.
- `docs/REPRODUCE.md` — step-by-step reproduction with the exact campaign grids, run counts
  and expected outcomes at every stage.
- `docs/DATA.md` — `.npz` naming convention, field-by-field archive schema, campaign
  inventory, versioning policy and known data caveats.
- `docs/CODE_MAP.md` — per-file catalogue of `code/`, classifying every module as live
  library, CLI entry point, manuscript tooling, optional-dependency tool, or historical
  scratch.
- `docs/LEGACY.md` — rationale for the dead and one-off scripts that are retained in the tree.
- `pack_repo.ps1` — builds a code-only upload archive and a full archive including raw data.

### Fixed

- **Portability:** `code/figures_ctp.py` hard-coded an absolute path to a personal tool
  directory (`C:\Users\<user>\.dsh\skills\nature-figure`). It now resolves through the
  `DSH_FIGURE_SKILL` environment variable with a home-relative default, and skips the audit
  gate cleanly when the helpers are absent. Figure rendering never depended on them.
- **Portability:** `code/figures_ctp.py` wrote preview rasters into a hard-coded
  `.dsh_tmp/figpreview` directory. The location is now the `BD_FIG_PREVIEW` environment
  variable, defaulting to `<repo>/.figpreview`. The stale `(+ .dsh_tmp preview)` console
  message was corrected to match.

### Changed

- Documentation now records that `run_campaign.py` rewrites `manifest_bd.json` on every
  invocation, so the shipped `spec` and `fgnspec` manifests are incomplete. The cell files,
  not the manifests, are documented as the authoritative grid inventory.

### Removed

- Nothing. No results, data, figures or analysed numbers were changed in this release. The
  changes are documentation, packaging and the two portability fixes above.

## [1.0.0] — 2026-10-09 — initial archived release

The state of the code as it produced the manuscript *Noise-generator bias in temporally
correlated ballistic-deposition simulations* (submitted to Physical Review E).

### Model and generators

- `code/bd_core.py`: (1+1)-dimensional ballistic deposition and RSOS rules under a parallel
  update with flat initial conditions and periodic boundaries; binarised temporally
  correlated noise.
- Three noise generators: Mandelbrot's FFGN (the published-literature generator, streaming
  and `O(N)` in memory) and two spectral filters — a bare power law `|f|^{-theta}` and the
  exact fractional-Gaussian spectral density.
- Counter-based splitmix64 RNG indexed by `(run, t, x)`, so results are deterministic and
  independent of thread count; numba `prange` parallel ensemble simulator.

### Generator bias

- Established that the FFGN normalised correlation amplitude is essentially
  `theta`-independent (`A_FFGN ≈ 1.27`) while the exact fGn amplitude `2 theta (theta + 1/2)`
  vanishes at `theta = 0`, so `FFGN(theta -> 0+)` is not white noise.
- Shown the effect is intrinsic by converging the generator in the number of AR components
  (`N = 23 / 32 / 48`).
- Binarised calibration against the exact fGn family, showing the entire nominal FFGN axis
  maps to `theta_eff >= 0.1888`, i.e. above the RG threshold `1/6`.

### Structural result

- Proved and verified `|Delta h| <= 1` from a flat initial condition, hence an increment pmf
  with support `{-1, 0, +1}`, exact `sum_x Delta h = 0` on a periodic lattice, `p_up = p_down`
  per configuration, and `skew(Delta h) ≡ 0` to machine precision.

### Diagnostics

- `code/diagnostics.py`: a null-calibrated Hellinger distance on **exact** discrete pmfs
  (binning-free), and an Isolation-Forest novelty detector on amplitude-free
  scaling/morphology features with a threshold calibrated at fixed false-positive rate.
- Head-to-head rejection, with data, of three of four candidate distributional observables.
- Onset extraction with bootstrap error bars, and onset finite-size-scaling fits both with a
  free `theta_c` and with `theta_c` forced to zero.

### Physics results

- The fixed-significance onset `theta*(L)` drifts downward over a factor of 16 in `L` as an
  artefact of the size-dependent split-half null floor; a forced-zero FSS fit is formally
  excellent, which demonstrates the construction cannot be trusted.
- The fixed-`theta` test: `DeltaH(theta, L)` is `L`-independent within errors for
  `theta <= 0.25` and significantly positive for `theta >= 0.30`, so the crossover lies
  beyond `L ≈ 512` for the smaller exponents. **No measured `theta_c` is reported.**
- Roughness exponents with run-bootstrap error bars, showing residual-based errors understate
  the uncertainty by roughly a factor of 8.
- A negative-`theta` control campaign (`theta = -0.2, -0.1, -0.05`) recovering `alpha = 1/2`.
- An exact-fGn-spectrum campaign closing the bare-power-law amplitude deficit as the
  systematic that could have faked a flat `alpha`; the two corrected filters agree to within
  0.008.
- Feature ablation and diagnostic-correlation analysis of the two diagnostics.

### Manuscript

- `manuscript/main.tex`: REVTeX 4.2, Physical Review E two-column layout; four full-width
  figures at exactly 178 mm; `manuscript/refs.bib`.
- `code/figures_ctp.py`: the four manuscript figures, exported as PDF (editable TrueType
  text) + SVG + 600 dpi TIFF, each passing a rendered-geometry collision audit.
- Static manuscript pre-flight in place of a local compile: `code/check_tex.py` and
  `code/table_width.py`.
- `code/build_overleaf.py`: Overleaf-ready build tree and upload archive.
