# Reproduce the 30 October 2026 PRE revision

The manuscript is `manuscript/main.tex`; the exported PDF is `manuscript/main.pdf`.
Use the matching source ZIP under `manuscript/` for Overleaf. It includes both
generated result-macro files and five PDF/SVG figures.

## Inputs

- `data/raw/rev20261026_shared`, `rev20261026_powerlaw`, `rev20261026_fgn`: 105 production cells, 64 independent runs each. Only these corrected campaigns are versioned. NPZ files contain the simulation specification, seed, widths, times and profiles. Spectral archives carry the corrected layout identifier.
- `data/analysis/revision_20261026/summary.json` and `noise_validation.json`: original revised statistics and field-validation results.
- `data/analysis/submission_20261030/paired.json`: 70 input archive hashes and seven paired slope estimates.
- `data/analysis/submission_20261030/summary.json`: paired inference and duration ratios, bootstrap settings, control-file hashes and protocol.
- Four `fixed_record_*.npz` files: three vectors of 32 per-run terminal mean squared widths, plus T, L, theta, seed and record length. They retain the estimator inputs; the binary fields can be regenerated from the recorded seed.

## Commands

Install Python scientific dependencies with `python -m pip install -r requirements.txt pytest`.

```sh
python -m pytest -q --basetemp .dsh_tmp/pytest-submission
python code/revision_analysis.py
python code/export_revision_tables.py
python code/revision_figures.py
python code/submission_controls.py
python code/build_overleaf.py
```

The supplemental script reuses its four cached control archives. To regenerate
the simulations independently, move those four NPZ files to a backup directory
and rerun it. They use 32 new independent runs at theta=0.10,0.40 and L=128,256.
Each run receives a length-4T fGn-spectrum record, whose nested T/2T/4T prefixes
are applied to flat initial interfaces. T=21870 and 61440 for L=128 and 256.
The terminal fifth of regularly sampled squared widths defines the estimator.
The low-frequency cutoff is fixed within each control cell; it differs from
the original production record's cutoff.

Paired Gaussian-filter bootstrap keeps run pairs together, draws independently
between sizes and uses 2000 replicates. The percentile intervals are pointwise,
not simultaneous confidence bands or prespecified equivalence tests.

Times New Roman regular, bold and italic fonts must be installed to regenerate
figures faithfully. All original four figures use `configure_style()` in
`revision_figures.py`; the supplemental figure uses the same configuration.
Compile the manuscript with REVTeX 4.2, BibTeX and a standard TeX distribution.

## Interpretation

The seven paired slope intervals include zero. All eight duration-ratio
intervals include one, and four-times-duration central ratios range from
0.989 to 1.034. These observations do not establish equivalence, saturation,
a critical temporal exponent or independence from cutoff/initialization.
Shared-versus-independent-innovation growth simulations were not added.

## Verification

All 50 tests passed with a project-local temporary directory. The new checks
compare nested-prefix simulations to a direct parallel-update reference and
verify supplemental central slopes against the original production fits.
An independent read-only audit reproduced the paired bootstrap and all archived
duration ratio summaries. It did not regenerate every input noise field.

The supplemental simulation environment used Python 3.12.14, NumPy 2.5.4,
SciPy 1.18.1, Numba 0.68.1, llvmlite 0.50.0 and Matplotlib 3.11.2.
Full regression tests were run with the existing project Python installation.
