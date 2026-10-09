# Legacy and historical scripts

`code/` contains scripts that are no longer part of the live pipeline. They are retained
deliberately, and this file records why — so that a future reader does not waste time trying
to run them, and does not "clean them up" and lose the provenance of the manuscript.

Nothing here affects any reported result. Every number and figure in the paper comes from the
live modules listed in [`CODE_MAP.md`](CODE_MAP.md) sections 1–4.

---

## 1. Dead scripts (reference only — they will fail)

These import functions that were removed from `bd_core.py` during development. They are kept
because they document the first exploratory stage of the project.

| Script | Broken import | Superseded by |
|---|---|---|
| `stage1_checks.py` | `generate_columns`, `_spectral_amplitude` | `validate_model.py`, `validate_fgn_spectrum.py` |
| `check_spectrum.py` | `_spectral_amplitude` | `validate_fgn_spectrum.py` |
| `debug_bd.py` | `_bd_ensemble` | `check_single_step.py`, direct `simulate()` calls |

`stage1_checks.py` is the original stage-1 scouting script: it validated the quantile-mapped
spectral filtering sequence, benchmarked one BD ensemble as a function of `L` and `T`, and did
the first KPZ sanity check at `theta = 0`. Its role was taken over by the dedicated validation
scripts, and its generator API was replaced when the noise generation moved to the
memory-bounded chunked implementation in `bd_core.generate_spectral_eta`.

## 2. One-off scratch utilities

These each performed a **single anchored edit** on a specific revision of
`manuscript/main.tex`. They exist to document exactly what was changed during the review and
formatting passes, and several contain the literal before/after strings. Re-running them
against the current manuscript would either fail an anchor-count assertion or corrupt the
file.

| Script | One-shot action |
|---|---|
| `fix_layout.py`, `fix_layout2.py` | Repaired overfull-box layout defects found by `overfull.py` in one specific PDF build. |
| `fix_tables.py` | Repaired two regressions introduced by the preceding layout fix. |
| `fix_inline_tables.py` | Rewrote seven specific inline tables to fit a single PRE column. |
| `transpose_tables.py` | Transposed two tables that remained too wide even after the rewrite. |
| `shorten_header.py` | Shortened one literal table header string. |
| `fix_author.py` | Removed author placeholder text, then rebuilt the Overleaf package. |
| `add_citations.py` | Added seven specific `refs.bib` entries that were defined but never cited. |
| `add_fgn_results.py` | Inserted the exact-spectrum supplementary paragraphs and Results subsection. |
| `finalize_fgn.py` | Updated the `fig:drift` caption for the three-generator version; contains the old caption text verbatim. |
| `_tmp_v7.py` | Cross-checked a v7 PDF against a formatting checklist; contains a hard-coded local path. |
| `_tmp_which_sample_count.py` | Answered one question about the effective sample count per cell. |

The `_tmp_*` prefix marks files that were scratch from the moment they were written. They are
kept only because they are referenced in the changelog history of the formatting passes.

**Release note.** If you are packaging this repository for public consumption and prefer a
tidier tree, moving `code/` sections 5 (one-off scratch) into a `code/legacy/` directory is
safe *as long as* you also check whether the script derives paths from `__file__`: most of
these do, and moving them one level deeper would make them read and write the wrong tree. If
in doubt, leave them in place — they cost nothing and they document real work.

## 3. Superseded but functional scripts

These still run correctly. They are not dead, just no longer used.

| Script | Status |
|---|---|
| `calibrate_generator.py` | First-generation calibration performed in the **Gaussian** domain. Superseded by `calibrate_binarised.py`, which calibrates after thresholding — the correct thing for this model, because thresholding compresses correlations via `C_bin = (2/pi) arcsin C`. Citing the Gaussian-domain `theta_eff` as the model's effective exponent is an error. |
| `figures.py`, `fig_generator.py`, `fig_exponents_calibrated.py`, `fig_onsets.py`, `fig_morphology.py`, `fig_drift.py` | Produced the earlier single-column draft figures now in `figs/generator/` and `figs/physics/`. Superseded by `figures_ctp.py` and `figs/ctp/`. Useful as a record of how the physics was first explored and plotted. |
| `zotero_probe.py` | A general Zotero library reader used to compile `manuscript/refs.bib`. Works, but needs `pymupdf` for the PDF-reading mode and defaults to `~/Zotero`. |

## 4. Not in this repository

The working tree also contained a scratch directory (`.dsh_tmp/`) with roughly two dozen
further one-off scripts — the v6/v7/v8 fix-and-verify pairs, figure-preview crops and
intermediate extraction dumps. Those were session working files, not part of the study, and
they are excluded from the repository by `.gitignore` and from the release archives by
`pack_repo.ps1`. Their effect on the manuscript is fully captured by the files listed in
section 2 above.
