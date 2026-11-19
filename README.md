# Temporally correlated ballistic deposition: data and code

This repository distributes numerical evidence and reproducibility software. Manuscript drafts, article PDFs, TeX/BibTeX sources and Overleaf bundles are not part of this data release.

## Numerical data

Fixed version: [data-20261119](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/tag/data-20261119).

- [Current numerical archive](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/download/data-20261119/bd_strict_binary_data_20261117.zip): parity controls, covariance perturbations, joint observation-time/record-period/mode controls, seeds, frozen producer sources and analysis inputs.
- [Historical innovation-control archive](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/download/data-20261119/innovation_controls_data_20261031.zip): calibration and validation records, including initial failed checks and refinements. It is also nested in the current archive.
- [SHA-256 checksums](https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/download/data-20261119/SHA256SUMS.txt).

The archives contain run-resolved NPZ arrays, JSON metadata and processed summaries, CSV figure-source data, simulation/analysis code, and reproduction instructions. They contain no article PDFs or TeX/BibTeX files. The numerical archive includes new simulations and historical inputs needed for the displayed comparisons; it does not redistribute every historical production record.

## Reproduce and inspect

Extract the current archive into an empty directory and start with [the numerical reproduction instructions](docs/strict_binary_20261117/REPRODUCE.md). Extract only missing inputs from the nested historical archive for the complete legacy audit; retain newer code and manifests.

- [Parity, residual-sensitivity and joint-limit summaries](data/analysis/strict_binary_20261117/summary.json).
- [Observation/period/mode control summaries](data/analysis/controlled_20261117/summary.json).
- Source CSVs and data-derived figures: `figs/strict_binary_20261117/` and `figs/refocus_20261117/`.
- Simulation, analysis and tests: `code/` and `tests/`.
- Parameters, seeds, resampling units and coverage limits: the reproduction instructions.

All 85 existing tests passed in the documented project Anaconda environment. Each numerical archive has an internal SHA-256 manifest; the release checksum file verifies the downloadable ZIP files.

## Cite the data

Yong Huang, *Numerical data and code for temporally correlated ballistic deposition*, version data-20261119, GitHub (2026), https://github.com/HuangYong2006/bd-correlated-noise-1.0.0-code/releases/tag/data-20261119.

Machine-readable citation metadata are in [CITATION.cff](CITATION.cff). This record has a fixed version URL; no DOI has been assigned. The existing repository license is in [LICENSE](LICENSE).

## Interpretation boundary

The parity comparison identifies a finite-scale response under exact binary two-time matching. Long-range covariance residuals and convergence limitations remain. These data do not establish a new universality class.
