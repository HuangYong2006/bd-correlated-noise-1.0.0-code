# Review notes index

These are the pre-submission review and formatting reports produced while preparing the
manuscript for Physical Review E. They are working documents, written in Chinese, and they
record the reasoning behind many of the formatting and presentation decisions in
`manuscript/main.tex`.

They are kept at the repository root rather than under `docs/` because code in the tree and
earlier review passes reference them by their root-relative filenames. This file indexes them
so the root directory does not need to be browsed to find them.

## Manuscript preparation

| File | Contents |
|---|---|
| [`PRE_review_edit_plan.md`](../PRE_review_edit_plan.md) | The edit plan derived from the first referee-style review: which claims to soften, which to support, and the order of the revisions. |
| [`PRE_review_current_status_v6.md`](../PRE_review_current_status_v6.md) | Status snapshot at revision v6, after the structural edits. |
| [`PRE_final_formatting_checklist_v7.md`](../PRE_final_formatting_checklist_v7.md) | Layout and formatting checklist for v7: float placement, table widths, page-break behaviour. |
| [`PRE_submission_readiness_v8.md`](../PRE_submission_readiness_v8.md) | Submission-readiness assessment at v8: page count, abstract length, metadata, and the list of blocking items (author metadata, empty `(Dated:)`, corresponding e-mail). |

**Status:** the blocking items listed in the v8 report were subsequently fixed in
`manuscript/main.tex` — the author is set, the hyperref metadata carries the real title and
author, and the placeholder date and acknowledgments text were removed. The remaining open
item is the corresponding e-mail address, which is a submission-system field.

## Physics and literature notes

| File | Contents |
|---|---|
| [`LR_temporal_correlation_numerics_report.md`](../LR_temporal_correlation_numerics_report.md) | Survey of the literature on long-range temporally correlated growth, and the numerics conventions used to reproduce it. |
| [`kpz_temporal_corr_theory_report.md`](../kpz_temporal_corr_theory_report.md) | Theory background: the RG threshold `theta_th = 1/6` and the expected exponent branches. |
| [`kpz_temporal_corr_exponent_table.md`](../kpz_temporal_corr_exponent_table.md) | Compilation of published exponent and `theta_c` values, reference by reference. |
| [`jiang_xia_2026_CTP_extraction.md`](../jiang_xia_2026_CTP_extraction.md) | Structured extraction of the most recent published study, used for the reference-by-reference implementation comparison. |

## Scratch and extraction dumps

`.dsh_tmp/` holds intermediate extraction dumps (literature PDFs converted to text, HTML
sources, and the v6/v7/v8 fix-and-verify scripts). It is working-session residue, excluded
from version control and from the release archives. Nothing in it is needed to reproduce any
result; see [`LEGACY.md`](LEGACY.md) section 4.
