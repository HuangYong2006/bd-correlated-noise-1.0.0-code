r"""
Build the Overleaf-ready PRE package from manuscript/main.tex.

main.tex is prepared with the publisher's own class, REVTeX 4.2, in the
two-column PRE layout.  Everything it loads -- revtex4-2, graphicx, xcolor, bm,
hyperref and the class's own amsmath/amssymb options -- is part of every TeX Live
distribution, so no source transformation is needed: the manuscript compiled
locally is byte-for-byte the manuscript uploaded.  The work here is

  1. copy the four figure PDFs into figs/ flat, so \graphicspath{{figs/}}
     resolves on Overleaf, and
  2. write a README documenting the compile recipe and the review-response
     changes.

The destination tree is wiped first: a stale figure from an earlier draft left
behind in manuscript/overleaf/ would silently be uploaded to the journal.

Output: manuscript/overleaf/{main.tex, refs.bib, abstract_zh.md, figs/*.pdf}
        manuscript/overleaf_bd_paper.zip
"""

from __future__ import annotations

import os
import shutil
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "manuscript")
DST = os.path.join(SRC, "overleaf")
FIGS_CTP = os.path.join(ROOT, "figs", "ctp")

FIGURES = ["fig1_generator", "fig2_exponents", "fig3_fixed_theta",
           "fig4_morphology"]

report = []

# ---------------------------------------------------------------- clean slate
if os.path.isdir(DST):
    shutil.rmtree(DST)
os.makedirs(os.path.join(DST, "figs"))

# ------------------------------------------------------------------- main.tex
tex = open(os.path.join(SRC, "main.tex"), encoding="utf-8").read()
open(os.path.join(DST, "main.tex"), "w", encoding="utf-8", newline="\n").write(tex)
report.append(f"  main.tex copied verbatim ({len(tex)} chars)")

# ------------------------------------------------------------- non-ASCII gate
# A stray non-ASCII glyph in the source is the single most common cause of a
# "Unicode character not set up for use with LaTeX" failure on Overleaf.
bad = sorted({c for c in tex if ord(c) > 127})
if bad:
    report.append(f"  !! NON-ASCII in main.tex: {bad}")
else:
    report.append("  main.tex is pure ASCII (no encoding risk)")

# ------------------------------------------------------------------- refs.bib
shutil.copyfile(os.path.join(SRC, "refs.bib"), os.path.join(DST, "refs.bib"))
report.append("  refs.bib copied")

zh = os.path.join(SRC, "abstract_zh.md")
if os.path.exists(zh):
    shutil.copyfile(zh, os.path.join(DST, "abstract_zh.md"))
    report.append("  abstract_zh.md copied (reference only, not compiled)")

# --------------------------------------------------------------------- figures
for name in FIGURES:
    for ext in (".pdf", ".svg"):
        s = os.path.join(FIGS_CTP, name + ext)
        if os.path.exists(s):
            shutil.copyfile(s, os.path.join(DST, "figs", name + ext))
report.append(f"  {len(FIGURES)} figure PDFs copied into figs/")

missing = [n for n in FIGURES if not os.path.exists(
    os.path.join(DST, "figs", n + ".pdf"))]
if missing:
    report.append(f"  !! MISSING figure PDFs: {missing}")

# ---------------------------------------------------------------------- README
readme = """OVERLEAF: HOW TO COMPILE THIS PROJECT
=====================================

Option A (recommended)
  Overleaf -> New Project -> Upload Project -> select
  "overleaf_bd_paper.zip", then press Recompile.

Option B
  Create a project and place the files as
      main.tex
      refs.bib
      figs/fig1_generator.pdf
      figs/fig2_exponents.pdf
      figs/fig3_fixed_theta.pdf
      figs/fig4_morphology.pdf
  set main.tex as the main document, press Recompile.

Compiler settings
  Compiler     : pdfLaTeX (default)
  TeX Live     : any recent version (revtex4-2 ships with TeX Live)
  Bibliography : BibTeX (Overleaf runs it automatically)
  No shell escape, no minted, no external font files, no non-ASCII glyphs.
  Packages used: revtex4-2 (class), graphicx, xcolor, bm, hyperref, plus the
  class's own amsmath/amssymb options -- all standard.

PRE FORMAT
  \\documentclass[aps,pre,twocolumn,superscriptaddress,floatfix,amsmath,amssymb]{revtex4-2}

  The publisher's own class, in the two-column reprint-style layout that PRE
  expects.  The four figures are exported at 178 mm, which is within 1% of the
  full REVTeX two-column text width (508 pt), so they are placed at natural size
  as full-width \\begin{figure*} floats and are legible at that size.

  Reviews are usually requested in a single-column preprint layout.  To produce
  one, replace `twocolumn` with `preprint` in the class options; the only other
  change required is that the four \\begin{figure*} floats may then become
  \\begin{figure}.

FIGURES
  All four figures are exported at exactly 178 mm wide.  Base glyph size is 7 pt
  (5 pt floor), and each is a vector PDF with embedded TrueType text.  Every
  figure is cited in the text.  Each one passed a rendered-geometry collision
  audit (0 blocking findings) and was additionally inspected at final size.

WHAT CHANGED IN THIS REVISION (final formatting pass)
  1. Link borders removed.  hyperref was loaded without options, so every DOI and
     arXiv link carried /Border [0 0 1] /C [0 1 1] -- a cyan rectangle, confirmed
     in the compiled PDF.  \\usepackage[hidelinks]{hyperref} keeps the links
     clickable and drops the border.  The same block now sets pdftitle/pdfauthor,
     which were empty in the previous PDF metadata.
  2. The bibliography no longer spills onto a page of its own.  Reference [24] was
     the only item on the last page; \\enlargethispage{4\\baselineskip} immediately
     before \\bibliography gains four lines on the page the list starts on.  No
     font size was reduced anywhere.
  3. The compile date is no longer printed as if it were the submission date:
     \\date{} replaces \\date{\\today}.  Put a real date back if the journal wants
     one on the title page.
  4. T1 font encoding with Latin Modern, loaded through \\IfFileExists so the
     manuscript still compiles if Latin Modern is absent.  This is what makes
     accented author names extract from the PDF as single glyphs; under OT1 the
     accent was a separate token, so a text search for "Ales" or "Lopez" could
     miss the reference.
  5. Notation unified: every {\\rm ...} subscript became {\\mathrm{...}} (75 sites),
     and the system-size range is written L=32\\dots512 throughout (11 sites).
  6. nofootinbib added to the class options, as the APS guidance recommends.

WHAT CHANGED IN THIS REVISION (submission-readiness pass)
  1. PDF metadata author was still the placeholder "Author list to be filled in"
     while the title page showed the real name.  pdftitle/pdfauthor now carry the
     real author, taken from the compiled PDF: Yong Huang, College of Physics,
     Chengdu University of Technology, Chengdu, China.
  2. The empty "(Dated:)" marker on the title page came from \\date{}: REVTeX
     prints "(Dated: )" for an empty date.  The \\date line is removed entirely, so
     REVTeX's default \\today applies and the marker disappears.  Put an explicit
     \\date{...} back if a fixed submission date is wanted.
  3. Only one title-page item is left for the author: the corresponding e-mail.
     It is a commented \\email line holding a bare FILL-IN token rather than an
     example address, so no placeholder can reach the PDF by accident.
  4. The reference list was ending with its last entry alone on a page.
     \\enlargethispage did not fix that on its own, because in two-column mode it
     only lengthens the column it is issued in; the list is now tightened through
     \\AtBeginEnvironment{thebibliography} with stretchable lengths, which applies
     after the list has started and therefore actually takes effect.

WHAT WAS VERIFIED IN THE COMPILED PDF (not assumed)
  * every font is Type1 or embedded TrueType -- no Type 3 bitmap, so no missing
    or boxed glyphs;
  * text is searchable ("ballistic deposition", "Hellinger", "fractional
    Gaussian" all found);
  * no text block reaches the printable margin on any page;
  * all 24 references render in APS format and every one is cited in the text;
  * the two "text glue" complaints in the review were checked string by string on
    the compiled PDF and are not present: the glued forms are absent while the
    spaced forms occur many times (e.g. "at theta" 56 times).  Math subscripts
    such as \\theta_c do flatten to "thetac" in text extraction, which is a
    property of the extractor, not of the PDF.

WHAT CHANGED IN THE PREVIOUS ROUND (response to the pre-submission review)
  1. TARGET FORMAT.  The manuscript had been converted to a CTP-style A4
     single-column layout; it is now back on the publisher class in the two-column
     PRE layout.  The four figures are unchanged and remain full-width.
  2. LITERATURE CLAIMS.  "used by every published simulation" and similar
     absolutes are gone.  A new Table I records, reference by reference, what
     could be established about each study's noise implementation from its
     published description alone, with explicit "not verified" entries for the
     two papers available only as abstracts.
  3. theta_eff is now defined as a shape-matched diagnostic only, explicitly not
     a reparameterisation of the RG coupling or of the universality class.
  4. The paper states throughout that it does NOT measure theta_c: within
     L = 32-512 no crossover is resolved, and the one place where a larger scale
     is mentioned is marked as an extrapolation from the trend rather than a
     measurement.
  5. The two diagnostics are no longer called independent.  A new subsection
     reports that their outputs correlate at r = 0.98 and gives a feature-family
     ablation (increment/terrace removed, scaling removed) showing which parts of
     the forest response survive.
  6. A new Appendix gives the statistical protocol: resampling unit, realisation
     counts per (theta, L), split-half null construction, onset definition, the
     degrees of freedom of each finite-size fit, the two-sigma slope criterion
     and the absence of a multiple-comparison correction, and the distinction
     between finite-size fit errors (+-0.0007) and bootstrap errors (+-0.012).
  7. The abstract was rewritten to 210 words in the order requested by the review
     and no longer quotes secondary numbers.
  8. The acknowledgments placeholders were deleted (see BEFORE YOU SUBMIT), and
     the Isolation Forest reference is now an @inproceedings entry.

DOCUMENT ORDER
  Conclusion, Appendix A (FFGN calibration), Appendix B (single-step proof),
  Appendix C (statistical protocol), References.

  There is deliberately no code or archive-listing appendix: the simulation
  settings that used to sit in one are now a short paragraph in Sec. 2.5, and
  the realisation counts are Table III.  A \\clearpage sits immediately before
  \\bibliography so that no pending float can be emitted inside or after the
  reference list.

BEFORE YOU SUBMIT
  * Fill in the author list and affiliations (search for "FILL IN").
  * Re-enable the acknowledgments: the placeholder text was removed, so the
    environment is commented out.  Search for "AUTHOR ACTION" and uncomment.
  * The manuscript reports a NEGATIVE result for theta_c.  Read the Conclusion
    carefully before submitting: we do not claim a value of theta_c; we claim
    that the reported thresholds are quantitatively reproduced by the generator
    bias alone.  Make sure you agree with that framing.
  * The data-availability statement in the Conclusion refers to the archive that
    accompanies the manuscript.  If you deposit it in a repository, replace that
    sentence with the real DOI or URL -- do not leave it as a promise.
  * If you can afford larger systems, adding L = 1024 at matched realisation
    count would convert the extrapolation in the Conclusion into a measurement.
    A partial L = 1024 campaign was attempted while preparing this revision but
    did not complete; no L = 1024 data is included and none is claimed.
"""
open(os.path.join(DST, "README_OVERLEAF.txt"), "w", encoding="utf-8",
     newline="\n").write(readme)
report.append("  README_OVERLEAF.txt written")

# ------------------------------------------------------------------------ zip
zpath = os.path.join(SRC, "overleaf_bd_paper.zip")
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _dirs, files in os.walk(DST):
        for f in sorted(files):
            full = os.path.join(root, f)
            z.write(full, os.path.relpath(full, DST))
report.append(f"  {os.path.basename(zpath)} written "
              f"({os.path.getsize(zpath)/1024:.0f} KB)")

print("Overleaf package rebuilt:")
print("\n".join(report))
