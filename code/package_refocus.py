"""Package verified paper sources and the completed, reproducible data revision."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT/'docs/research_refocus_20261117'
STRICT = ROOT/'docs/strict_binary_20261117'
PAPER = ROOT/'manuscript'
DIST = ROOT/'dist'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_package(target, files, readme):
    manifest = []
    with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, source in sorted(files.items()):
            data = source.read_bytes()
            archive.writestr(name, data)
            manifest.append(dict(file=name, bytes=len(data), sha256=digest(data)))
        archive.writestr('README_PACKAGE.txt', readme)
        archive.writestr('MANIFEST.json', json.dumps(manifest, indent=2)+'\n')
    with zipfile.ZipFile(target) as archive:
        assert archive.testzip() is None
        for row in manifest:
            assert digest(archive.read(row['file'])) == row['sha256']
    return dict(file=target.relative_to(ROOT).as_posix(), members=len(manifest)+2,
                bytes=target.stat().st_size, sha256=digest(target.read_bytes()))


def package_paper(target):
    names = [
        'main.tex', 'supplemental.tex', 'refs.bib',
        'figs/fig1_noise_statistics.pdf', 'figs/fig2_widths.pdf',
        'figs/fig3_effective_flow.pdf', 'figs/fig4_strict_binary.pdf',
        'figs/fig5_convergence.pdf', 'figs/figS1_noise_validation.pdf',
        'figs/fig5_controls.pdf', 'figs/fig3_fixed_theta.pdf',
        'figs/fig4_morphology.pdf',
    ]
    files = {name: PAPER/name for name in names}
    return write_package(target, files,
        'Physical Review E submission sources. Compile main.tex and supplemental.tex '
        'independently with pdfLaTeX/BibTeX (REVTeX 4.2) or Tectonic. Each document '
        'has a single TeX file. Upload the compiled article and supplemental PDFs separately. '
        'MANIFEST.json verifies the source files. '
        'The repository stores the numerical archive in dist/bd_strict_binary_data_20261117.zip.\n')


def main():
    audit = json.loads((DOCS/'campaign_audit.json').read_text())
    assert audit['passed'] and audit['independent_interface_trajectories'] == 640
    summary = json.loads((ROOT/'data/analysis/controlled_20261117/summary.json').read_text())
    assert summary['extension']['runs_per_size'] == [64]*5
    paper = package_paper(PAPER/'overleaf_bd_strict_binary_20261117.zip')
    data = {}
    def add(path):
        data[path.relative_to(ROOT).as_posix()] = path
    for group in ['baseline', 'period4', 'period8', 'modes_plus4', 'calibration_sensitivity']:
        for path in (ROOT/'data/raw/controlled_20261117'/group).rglob('*'):
            if path.suffix in {'.npz', '.json'}:
                add(path)
    strict_summary = ROOT/'data/analysis/strict_binary_20261117/summary.json'
    for rel, sha in json.loads(strict_summary.read_text())['sources'].items():
        assert digest((ROOT/rel).read_bytes()) == sha
    assert len(json.loads(strict_summary.read_text())['joint'])==6
    for path in (ROOT/'data/raw/strict_binary_20261117').rglob('*'):
        if path.suffix in {'.npz', '.json'}:
            add(path)
    for pattern in ['code/*.py', 'figs/strict_binary_20261117/*',
                    'data/analysis/strict_binary_20261117/*.json',
                    'docs/superpowers/plans/2026-11-17-strict-binary-control.md', 'tests/*.py', 'figs/refocus_20261117/*',
                    'figs/revision_20261026/*',
                    'data/analysis/controlled_20261117/*.json']:
        for path in ROOT.glob(pattern):
            if path.is_file():
                add(path)
    for source in json.loads((DOCS/'figure_provenance.json').read_text())['sources']:
        path = ROOT/source['file']
        assert digest(path.read_bytes()) == source['sha256']
        add(path)
    for path in DOCS.iterdir():
        if path.suffix in {'.md', '.json', '.py', '.xml'} and path.name != 'release_manifest.json':
            add(path)
    for path in STRICT.iterdir():
        if path.suffix in {'.md','.json','.xml'} and path.name != 'release_manifest.json':
            add(path)
    # Preserve the small historical assets required by the existing regression suite.
    for manifest in ['figure_sources.json', 'fig5_sources.json']:
        record = json.loads((ROOT/'figs/revision_20261026'/manifest).read_text())
        for key in ['summary', 'noise_validation', 'paired_sources']:
            if key in record:
                add(ROOT/record[key])
    historical = ROOT/'data/analysis/submission_20261030/summary.json'
    add(historical)
    for row in json.loads(historical.read_text())['duration']:
        add(ROOT/row['raw_file'])
    add(DOCS/'reproduced_analysis/summary.json')
    add(ROOT/'docs/revision_20261026/noise_validation_notes.md')
    for name in ['pyproject.toml', 'requirements.txt', 'LICENSE']:
        if (ROOT/name).is_file():
            add(ROOT/name)
    add(DIST/'innovation_controls_data_20261031.zip')
    numerical = write_package(DIST/'bd_strict_binary_data_20261117.zip', data,
        'Strict binary control revision. Start with docs/strict_binary_20261117/REPRODUCE.md. '
        'Extract only missing members of nested dist/innovation_controls_data_20261031.zip into this root when '
        're-running the complete historical stationary audit; it retains initial failed checks. '
        'The new campaign and all main-figure inputs are included directly. '
        'Not all 105 original production raw cells are redistributed; their frozen summaries are retained. '
        'The calibration refits share 16 innovations and are not 1040 independent records. '
        'The timing-only baseline_fast directory is excluded from inference and this package.\n')
    (STRICT/'release_manifest.json').write_text(json.dumps(dict(packages=[paper, numerical]), indent=2)+'\n')
    print(json.dumps([paper, numerical], indent=2))


if __name__ == '__main__':
    main()
