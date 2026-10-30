"""Build the revised PRE source package from an explicit file list.

Only listed sources are zipped, so stale files cannot enter the new package.
Original raw data and the separate pre-revision archive are preserved.
"""
from pathlib import Path
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'manuscript'
DST = SRC / 'overleaf'
FIGS = ROOT / 'figs' / 'revision_20261026'
FIGURES = ['fig1_generator', 'fig2_exponents', 'fig3_fixed_theta', 'fig4_morphology', 'fig5_controls']

README = '''OVERLEAF SOURCE PACKAGE - PRE REVISION 2026-10-30
=====================================================
Upload overleaf_bd_paper.zip as an Overleaf project and compile main.tex using
pdfLaTeX and BibTeX. REVTeX 4.2 and the packages used are standard TeX Live
components. revision_results.tex and submission_results.tex supply numerical tables
and must stay beside main.tex and refs.bib.

This package uses a regenerated 105-cell comparison: three inputs, seven theta
values, five sizes, 64 independent runs per cell. Five figures use corrected spectral data, including paired-filter
inference and fixed-record duration controls. All figure text is Times New Roman. The text distinguishes finite-size late-time slopes from
asymptotic exponents. It neither reports theta_c nor attributes errors to
previous authors. The former flat corrected KPZ-branch claim is withdrawn.

Files
  main.tex, refs.bib, revision_results.tex, submission_results.tex
  abstract_zh.md (reading aid, not part of the compiled manuscript)
  figs/fig1_generator through fig5_controls, each PDF and SVG

The standalone overleaf_bd_paper__1_.pdf elsewhere in the project is an old
draft, not a rendering of this package. Recompile these sources.
Code/data reproduction is described in docs/REPRODUCE.md in the project.
Before actual submission the author must verify the revised interpretation
and supply the real corresponding email. Compilation alone does not certify
scientific readiness or acceptance.
'''


def main():
    DST.mkdir(exist_ok=True)
    (DST / 'figs').mkdir(exist_ok=True)
    files = {name: SRC / name for name in
             ['main.tex', 'refs.bib', 'revision_results.tex', 'submission_results.tex', 'abstract_zh.md']}
    for name in FIGURES:
        for extension in ['pdf', 'svg']:
            files[f'figs/{name}.{extension}'] = FIGS / f'{name}.{extension}'
    for name, source in files.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, DST / name)
    (DST / 'README_OVERLEAF.txt').write_text(README, encoding='utf-8')
    zip_path = SRC / 'overleaf_bd_paper.zip'
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in sorted([*files, 'README_OVERLEAF.txt']):
            archive.write(DST / name, arcname=name)
    print(f'Rebuilt {zip_path} ({zip_path.stat().st_size / 1024:.0f} KB)')
    print(f'Included {len(files)+1} explicitly verified source files')


if __name__ == '__main__':
    main()
