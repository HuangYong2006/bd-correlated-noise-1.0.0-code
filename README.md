# Temporally correlated ballistic deposition

This project studies parallel-update ballistic deposition with explicitly specified
shared-innovation AR and Gaussian spectral noise inputs. The revised manuscript is
**Noise construction and finite-size inference in temporally correlated ballistic deposition**,
prepared in Physical Review E format.

## PRE revision, 30 October 2026

The current manuscript and compiled PDF are `manuscript/main.tex` and
`manuscript/main.pdf`. The five figures use Times New Roman throughout.
The matching Overleaf source archive is `manuscript/overleaf_bd_paper.zip`.

This revision adds a paired bootstrap for the two common-stream Gaussian
filters and four fixed-record duration-control cells (32 independent runs each).
All seven pointwise slope-difference intervals include zero; all eight duration
ratio intervals include one. These are non-detections at the available precision,
not evidence of equivalence, saturation, or an asymptotic critical parameter.

All 105 corrected production cells are versioned under the three
`data/raw/rev20261026_*` directories. Other raw campaigns remain excluded.
See [the current reproduction guide](docs/submission_20261030/REPRODUCE.md)
and [revision notes](docs/submission_20261030/修改说明.md).
The corresponding-author email is intentionally left for the author to fill.

## Scientific correction, 26 October 2026

A time/site array-layout error was found in the former production spectral path.
It is fixed and covered by end-to-end temporal/spatial correlation tests. The old
`spec` and `fgnspec` results are available only in Git history and are **not valid
corrected-noise evidence**. Earlier readiness notes and generator-bias conclusions
are superseded by the revision described here.

New tracks `rev20261026_shared`, `rev20261026_powerlaw` and `rev20261026_fgn` each
contain 35 cells: seven theta values, five sizes (32 through 512), and 64 independent
runs per cell. Their 105 cells were regenerated rather than relabelled. Spectral
archives carry `noise_layout=spectral-run-time-site-v2`; new analysis rejects older
spectral layouts.

The new data show increasing finite-size late-time width slopes and decreasing
flat-pair fractions for all three inputs. At theta=0.25, shared AR gives
0.717 +/- 0.008 and the fGn-spectrum filter 0.678 +/- 0.007. These are finite-window
slopes with run-bootstrap uncertainty, not measured asymptotic exponents.
No theta_c is reported and no historical implementation is identified as erroneous.
The old claims of a flat corrected KPZ branch, opposite-sign morphology and a
quantitative explanation of published thresholds have been withdrawn.

## Reproduction and files

Use [docs/REPRODUCE.md](docs/REPRODUCE.md) for current commands and
[docs/DATA.md](docs/DATA.md) for archive fields and legacy-data boundaries.

- `code/bd_core.py`: simulation and corrected space-time spectral assembly.
- `code/validate_production_noise.py`: final-interface correlation validation.
- `code/revision_analysis.py`: run-level widths, effective slopes and local-pmf distances.
- `code/export_revision_tables.py`: data-derived LaTeX numbers.
- `code/revision_figures.py`: the original four revised figures and source manifest.
- `code/submission_controls.py`: paired filter inference, duration controls and the fifth figure.
- `manuscript/main.tex`: current manuscript; `revision_results.tex` and `submission_results.tex` are required inputs.
- `manuscript/overleaf_bd_paper.zip`: rebuilt upload package.
- `docs/revision_20261026/`: correction, test and figure QA records.

```powershell
& 'D:/TOOLS/anaconda3/python.exe' -m pytest tests -q --basetemp .dsh_tmp/pytest
& 'D:/TOOLS/anaconda3/python.exe' code/validate_production_noise.py
& 'D:/TOOLS/anaconda3/python.exe' code/revision_analysis.py
& 'D:/TOOLS/anaconda3/python.exe' code/export_revision_tables.py
& 'D:/TOOLS/anaconda3/python.exe' code/revision_figures.py
& 'D:/TOOLS/anaconda3/python.exe' code/submission_controls.py
& 'D:/TOOLS/anaconda3/python.exe' code/build_overleaf.py
```

Obsolete scripts, old analysis outputs and duplicate source folders were removed
from the current tree. See [cleanup notes](docs/CLEANUP.md); historical versions
remain recoverable from Git history.
The project remains MIT-licensed. It does not include third-party full texts in
release packages. Publication metadata and data-repository identifiers must be
updated only when an actual submission or deposit exists.
