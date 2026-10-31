"""Self-contained numerical supplement and a manuscript/source/data bundle."""
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def build_data_archive(summary, result, sizes):
    payload = {'summary.json': summary.read_bytes()}
    for source in result['sources']:
        p = ROOT/source['file']
        data = p.read_bytes()
        if hashlib.sha256(data).hexdigest() != source['sha256']:
            raise ValueError(f'changed numerical source: {p}')
        payload[source['file']] = data
    # Retain baseline controls, including both original failures, for audit
    # and for regenerating the documented refinement.
    for size in sizes:
        for theta in [.1, .25, .4]:
            p = ROOT/'data/raw/innovation_binary_20261031_v3'/f'theta{theta:.2f}_L{size}.npz'
            payload[p.relative_to(ROOT).as_posix()] = p.read_bytes()
    scripts = ['bd_core', 'innovation_controls', 'innovation_streaming', 'control_archives',
               'binary_covariance_control', 'innovation_analysis', 'refine_binary_controls',
               'verify_submission_data', 'validate_large_binary_control',
               'update_innovation_manuscript', 'submission_archive', 'revision_figures',
               'independent_supplement_audit']
    for name in scripts:
        payload[f'code/{name}.py'] = (ROOT/'code'/f'{name}.py').read_bytes()
    for p in (ROOT/'docs/innovation_20261031/source_snapshot').glob('*'):
        payload['source_snapshot/'+p.name] = p.read_bytes()
    for name in ['requirements.txt', 'LICENSE']:
        payload[name] = (ROOT/name).read_bytes()
    payload['refinement_report.json'] = (ROOT/'docs/innovation_20261031/refinement_report.json').read_bytes()
    if 512 in sizes:
        report = ROOT/'docs/innovation_20261031/large_validation_report.json'
        if report.exists():
            payload['additional_validation_report.json'] = report.read_bytes()
    seeds = dict(source_cases=[], baseline_binary=[], recalibration=[],
                 additional_covariance_only=[],
                 bootstrap=dict(source_slopes=20261031, covariance_validation=20261031,
                                heldout_uniform_slopes=20261032, surrogate_slopes=20261033))
    for size in sizes:
        for theta in [.1, .25, .4]:
            name = f'theta{theta:.2f}_L{size}.npz'
            with np.load(ROOT/'data/raw/innovation_20261031_v3'/name, allow_pickle=False) as z:
                spec = json.loads(str(z['spec']))
            source_pairs, binary_pairs = [], []
            for run in range(spec['runs']):
                children = np.random.SeedSequence([spec['seed'], int(round(theta*1e6)), size, run]).spawn(2)
                source_pairs.append([int(c.generate_state(1)[0]) for c in children])
                children = np.random.SeedSequence([20261101, int(round(theta*1e6)), size, run]).spawn(2)
                binary_pairs.append([int(c.generate_state(1)[0]) for c in children])
            seeds['source_cases'].append(dict(theta=theta, L=size, root_seed=spec['seed'],
                recipe='SeedSequence([root_seed, round(theta*1e6), L, run]).spawn(2); each child.generate_state(1)[0]',
                gaussian_seed32=[p[0] for p in source_pairs], uniform_seed32=[p[1] for p in source_pairs],
                uniform_independent_counter_base='uint64(uniform_seed32) * 0x9E3779B97F4A7C15'))
            seeds['baseline_binary'].append(dict(theta=theta, L=size, root_seed=20261101,
                recipe='SeedSequence([20261101, round(theta*1e6), L, run]).spawn(2); each child.generate_state(1)[0]',
                gaussian_seed32=[p[0] for p in binary_pairs], uniform_reference_seed32=[p[1] for p in binary_pairs]))
    for size in [32, 64]:
        seeds['recalibration'].append(dict(theta=.4, L=size, root_seed=20261103,
            recipe='SeedSequence([20261103, L, domain, run]).generate_state(1)[0]; domains0/1/2=training/Gaussian/uniform reference',
            training_uniform_seed32=[int(np.random.SeedSequence([20261103,size,0,r]).generate_state(1)[0]) for r in range(128)],
            gaussian_validation_seed32=[int(np.random.SeedSequence([20261103,size,1,r]).generate_state(1)[0]) for r in range(64)],
            uniform_validation_seed32=[int(np.random.SeedSequence([20261103,size,2,r]).generate_state(1)[0]) for r in range(64)]))
    if 512 in sizes:
        seeds['additional_covariance_only'].append(dict(theta=.4, L=512, root_seed=20261105,
            recipe='SeedSequence([20261105,512,run]).spawn(2); default_rng(child0), child1.generate_state(1)[0]',
            gaussian_seed_sequences=[dict(entropy=[20261105,512,r], spawn_key=[0]) for r in range(128)],
            uniform_seed32=[int(np.random.SeedSequence([20261105,512,r]).spawn(2)[1].generate_state(1)[0]) for r in range(128)]))
    payload['SEEDS.json'] = json.dumps(seeds, indent=2).encode('utf-8')
    refinements = json.loads(payload['refinement_report.json'])
    index = dict(
        raw_npz_folders={
            'data/raw/innovation_20261031_v3':'15 six-variant source ensembles; spec, late_w2, profiles and input statistics',
            'data/raw/innovation_binary_20261031_v3':'15 original surrogate archives, including all original failed validations',
            'data/raw/innovation_binary_refined_20261031':'15 final surrogate archives; two recalibrations and the separate L512 validation'},
        original_failed_validations=[dict(theta=r['spec']['theta'], L=r['spec']['L'],
            file=r['previous_archive'], upper95=r['previous_gate']['simultaneous_upper95'],
            arrays=['sign_covariance','validation_reference_covariance','validation_lags'])
            for r in refinements],
        recalibrations=[dict(theta=r['spec']['theta'], L=r['spec']['L'],
            file='data/raw/innovation_binary_refined_20261031/'+Path(r['previous_archive']).name,
            arrays=['training_sign_covariance','fitted_spectrum','sign_covariance','validation_reference_covariance'],
            metadata='metadata JSON contains fit coefficients, source hashes, sampling protocol and root seed')
            for r in refinements],
        seed_list='SEEDS.json', independent_check='python code/independent_supplement_audit.py',
        analysis='code/innovation_analysis.py', figure='code/update_innovation_manuscript.py --summary summary.json --figure-only')
    if 'additional_validation_report.json' in payload:
        extra = json.loads(payload['additional_validation_report.json'])
        index['original_failed_validations'].append(dict(theta=.4,L=512,file=extra['previous_archive'],
            upper95=extra['previous_gate']['simultaneous_upper95'],
            arrays=['sign_covariance','validation_reference_covariance','validation_lags']))
        index['additional_covariance_only'] = dict(file='data/raw/innovation_binary_refined_20261031/theta0.40_L512.npz',
            arrays=['independent_validation_gaussian_covariance','independent_validation_uniform_covariance'],
            fitted_spectrum_and_interfaces='unchanged from the original archive',
            protocol='additional_validation_report.json')
    payload['DATA_INDEX.json'] = json.dumps(index, indent=2).encode('utf-8')
    sizes_arg = ' '.join(map(str, sizes))
    binary_folder = Path(next(s['file'] for s in result['sources']
                              if 'innovation_binary_refined_' in s['file'])).parent.as_posix()
    unchanged = sum(
        (ROOT/s['file']).read_bytes() == (ROOT/'data/raw/innovation_binary_20261031_v3'/Path(s['file']).name).read_bytes()
        for s in result['sources'] if 'innovation_binary_refined_' in s['file'])
    readme = f"""NUMERICAL SUPPLEMENT - INNOVATION CONTROLS

Completed sizes: {sizes_arg}; theta=0.10,0.25,0.40.
Each interface ensemble has 64 independent runs; binary-surrogate comparisons
use 32 held-out uniform-shared source runs. All seed/period parameters and
calibration metadata are in the NPZ archives. SEEDS.json explicitly lists
per-run seeds and their derivations; DATA_INDEX.json locates the original
failures, recalibrations and exact data fields. No pickle objects are needed.

Unzip this archive into an empty directory, then:
  python -m pip install -r requirements.txt
  python code/verify_submission_data.py
  python code/independent_supplement_audit.py
  python code/innovation_analysis.py --source data/raw/innovation_20261031_v3 --binary-source {binary_folder} --sizes {sizes_arg} --out reproduced_analysis
  python code/update_innovation_manuscript.py --summary summary.json --figure-only

The final command recreates paper Figure7 as PDF/SVG/PNG without a TeX engine.
It preserves the numerical content; local font availability may affect typography.

The verifier checks every payload hash, recomputes slope/bootstrap contrasts
and independent covariance-validation bounds from the raw run arrays, and
reports which cells pass the unchanged0.005 threshold.
The independent audit imports no project generator or analysis functions.
It recomputes the original failures, current bounds, least-squares slopes,
bootstrap intervals and four-sign moments, reconstructs recalibrated spectra,
and regenerates a calibration record from its explicit seed for each refined
cell. It writes INDEPENDENT_AUDIT.json for inspection.

To regenerate the two refined cells from independent noise, run:
  python code/refine_binary_controls.py
This uses 128 fresh 64-site calibration records and 64 new 64-site validation
records per failed cell, fixed seed root 20261103 and separate seed domains.
The fit includes every lag 0..64 plus the original long lags. There are no
further fit retries using these fresh validation results and no tolerance relaxation.

The original failed archives remain under data/raw/innovation_binary_20261031_v3.
The refined archives and their metadata identify the original file and hash.
There are {unchanged} byte-identical baseline archives in the selected grid.
The calibration records and fresh uniform covariance-reference records were
not used to fit interface exponents. Gaussian surrogate records supply both
interface trajectories and their baseline Gaussian covariance diagnostics.
Any additional covariance-only cohort is separate and supplies no interface
trajectories; its fixed-spectrum protocol and original failed archive are
identified in the cell metadata.

For the full five-size supplement, the additional L512 covariance-only cohort
can be regenerated with:
  python code/validate_large_binary_control.py
Its fixed fitted spectrum and original interface arrays are retained; 128 new
64-site Gaussian records and 128 independent 64-site uniform records verify
the same 0.005 criterion without changing fitted parameters or interface data.

Source snapshots record the exact earlier numerical producers. The current
scripts preserve their numerical model and add atomic file publication.
To regenerate a source campaign rather than analyze the supplied archives,
use a fresh output directory with code/innovation_controls.py; never overwrite
the archived evidence. Exact hashes identify archived bytes, while numerical
reproduction on other library versions should be assessed with tolerances.

summary.json contains 2000 whole-run bootstrap intervals. Inference is
conditional on the fitted surrogate spectrum; finite-lag covariance matching
does not establish equality of complete processes or asymptotic exponents.
"""
    payload['README_REPRODUCE.txt'] = readme.encode('utf-8')
    deposit_readme = (
        'Title: Noise construction and finite-size inference in temporally correlated ballistic deposition\n'
        'Authors: Yong Huang\n'
        'Description: Raw run archives and reproducible code for shared/independent innovations, '
        'covariance compensation, independent binary-covariance validation and four-time noise statistics.\n'
        'Numerical archive: innovation_controls_data_20261031.zip\n'
        'Contents: baseline and refined archives, original validation failures, calibration and validation metadata, '
        'seeds, exact producer snapshots, analysis and Figure 7 code, file hashes and standalone verifier.\n'
        'Requirements: Python 3 and the packages in requirements.txt.\n'
        'Instructions: Unzip into an empty directory and follow README_REPRODUCE.txt.\n')
    payload['README.TXT'] = deposit_readme.encode('utf-8')
    (ROOT/'manuscript/README.TXT').write_text(deposit_readme, encoding='utf-8')
    manifest = dict(format='innovation-supplement-v1', sizes=sizes,
                    entries=[dict(file=name, sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
                             for name, data in sorted(payload.items())])
    payload['MANIFEST.json'] = json.dumps(manifest, indent=2).encode('utf-8')
    path = ROOT/'manuscript/innovation_controls_data_20261031.zip'
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(payload.items()):
            archive.writestr(name, data)
    return path


def build_submission_bundle():
    names = ['BD_noise_controls_updated.pdf', 'overleaf_bd_paper.zip',
             'innovation_controls_data_20261031.zip', 'README.TXT']
    path = ROOT/'manuscript/PRE_submission_materials_20261031.zip'
    records = []
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in names:
            p = ROOT/'manuscript'/name
            records.append(dict(file=name, bytes=p.stat().st_size,
                                sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
            archive.write(p, name)
        archive.writestr('MANIFEST.json', json.dumps(records, indent=2))
        archive.writestr('README_SUBMISSION.txt',
                         'Manuscript PDF, editable LaTeX source package, and numerical supplement.\n'
                         'Provide innovation_controls_data_20261031.zip with the manuscript as '
                         'supplementary data/code. It includes original failures, refinement, '
                         'raw run arrays, code, seeds, hashes and a standalone verifier.\n'
                         'The completed control-size scope is stated in the manuscript and supplement.\n')
    return path
