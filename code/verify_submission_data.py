"""Verify an extracted numerical supplement without the author's workspace."""
import hashlib
import json
from pathlib import Path

import numpy as np

from innovation_analysis import bootstrap_slopes, covariance_gate, four_sign_moment, covariance_validation_samples
from innovation_controls import VARIANTS


def verify(root):
    manifest = json.loads((root/'MANIFEST.json').read_text(encoding='utf-8'))
    for entry in manifest['entries']:
        p = root/entry['file']
        assert hashlib.sha256(p.read_bytes()).hexdigest() == entry['sha256'], entry['file']
    summary = json.loads((root/'summary.json').read_text(encoding='utf-8'))
    archives = {}
    for entry in summary['sources']:
        p = root/entry['file']
        assert hashlib.sha256(p.read_bytes()).hexdigest() == entry['sha256']
        with np.load(p, allow_pickle=False) as z:
            if 'spec' in z.files:
                spec = json.loads(str(z['spec']))
                archives.setdefault((spec['theta'], spec['L']), {})['source'] = {
                    k: z[k] for k in z.files}
            else:
                metadata = json.loads(str(z['metadata']))
                spec = metadata['spec']
                archives.setdefault((spec['theta'], spec['L']), {})['binary'] = {
                    k: z[k] for k in z.files}
    passed = 0
    for cell in summary['cells']:
        pair = archives[cell['theta'], cell['L']]
        b = pair['binary']
        metadata = json.loads(str(b['metadata']))
        histories = [metadata, metadata.get('additional_covariance_validation', {})]
        for history in histories:
            if 'previous_archive' in history:
                old = root/history['previous_archive']
                assert hashlib.sha256(old.read_bytes()).hexdigest() == history['previous_archive_sha256']
        if 'additional_covariance_validation' in metadata:
            with np.load(root/metadata['additional_covariance_validation']['previous_archive'],
                         allow_pickle=False) as old:
                for key in ['fitted_spectrum', 'late_w2', 'profiles']:
                    np.testing.assert_array_equal(b[key], old[key])
        gate = covariance_gate(*covariance_validation_samples(b))
        assert gate['margin'] == .005
        np.testing.assert_allclose(gate['simultaneous_upper95'],
                                   cell['binary_covariance_gate']['simultaneous_upper95'], atol=1e-14)
        passed += gate['passed']
        a = four_sign_moment(b['persistence'], b['sign_covariance'], b['validation_lags'])
        ref = four_sign_moment(b['validation_reference_persistence'],
                               b['validation_reference_covariance'], b['validation_lags'])
        np.testing.assert_allclose(a.mean()-ref.mean(), cell['four_sign_moment_difference'], atol=1e-12)
    sizes = manifest['sizes']
    contrasts = 0
    for theta in [.1, .25, .4]:
        inputs = [archives[theta, size]['source']['late_w2'] for size in sizes]
        central, boot = bootstrap_slopes(sizes, inputs)
        for i, variant in enumerate(VARIANTS):
            row = next(r for r in summary['slopes'] if r['theta']==theta and r['variant']==variant)
            np.testing.assert_allclose(central[i], row['alpha'], atol=1e-12)
            np.testing.assert_allclose(np.quantile(boot[:, i], [.025, .975], axis=0), row['ci95'], atol=1e-12)
        ref = [archives[theta, s]['source']['late_w2'][32:, 3] for s in sizes]
        target = [archives[theta, s]['binary']['late_w2'] for s in sizes]
        a, ab = bootstrap_slopes(sizes, ref, seed=20261032)
        b, bb = bootstrap_slopes(sizes, target, seed=20261033)
        row = next(r for r in summary['comparisons'] if r['theta']==theta
                   and r['comparison']=='binary_match_vs_heldout_uniform')
        np.testing.assert_allclose(b-a, row['delta_alpha'], atol=1e-12)
        np.testing.assert_allclose(np.quantile(bb-ab, [.025, .975], axis=0), row['ci95'], atol=1e-12)
        contrasts += 1
    return dict(files_verified=len(manifest['entries']), cells_verified=len(archives),
                covariance_passed=int(passed), covariance_tolerance=.005,
                gaussian_uniform_slope_ensembles=18, surrogate_contrasts_verified=contrasts)


if __name__ == '__main__':
    print(json.dumps(verify(Path(__file__).resolve().parents[1]), indent=2))
