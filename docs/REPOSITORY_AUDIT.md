# Repository completeness audit, 30 October 2026

## Evidence completeness

- All 105 corrected production NPZ cells and four fixed-record control NPZ files are present and versioned.
- All 183 input-file/hash references in the production, paired and duration summaries resolved and matched. This count includes repeated references; it is not a count of distinct data files.
- The current five figures have all 15 PDF/SVG/PNG outputs. The two figure-source manifests now cover them with repository-relative paths and input hashes.
- The 16-entry Overleaf ZIP passed its CRC check. Its manuscript, bibliography, numerical macros, reading abstract and ten PDF/SVG figures match the canonical sources after normalizing text line endings; binary figures match byte for byte.
- No retained Python code or regression test imports a removed local module. The main reproduction sequence and the source-package builder are present.
- The canonical manuscript source and exported PDF were not changed by this audit.

## Gaps corrected

1. Added `submission_controls.py` to the main README and production reproduction command sequence. The detailed supplemental guide already contained it.
2. Replaced machine-specific absolute paths in the four-figure source manifest with repository-relative paths for repository inputs/outputs.
3. Added `fig5_sources.json` for the supplemental figure and a regression test covering all five figures' provenance.
4. Removed 11 further files that belonged to retired data/shape-calibration workflows, obsolete physical-threshold interpretations or personal literature tooling. The exact list is in `CLEANUP.md`; earlier versions remain in Git history.
5. Normalized figure-manifest text hashes and source-package text to UTF-8/LF so Windows and Unix checkouts agree. Binary archive hashes remain byte-exact.
6. Updated retained exploratory-tool descriptions and documentation to avoid references to removed tools or superseded interpretations.

## Intentionally excluded local material

Untracked pre-review folders, pytest temporary directories, compilation diagnostics,
the duplicate working manuscript body and the old source-archive directory are
working/history material, not missing inputs to the current manuscript. Ignored
temporary dependencies, caches, third-party literature and older raw campaigns
are likewise not publication evidence. They were not uploaded or deleted during
this audit.

The remaining optional PDF inspection, benchmark, model/field validation and
archive-integrity tools are usable supporting utilities, rather than missing
parts of the primary inference pipeline. `features.py`, `analyze.py`,
`diagnostics.py` and `feature_ablation.py` remain for tested run-level diagnostics.

## Verification

All 51 regression tests passed with a repository-local temporary directory.
The audit checks file completeness, source dependencies and provenance; it does
not establish scientific novelty or guarantee journal acceptance.
