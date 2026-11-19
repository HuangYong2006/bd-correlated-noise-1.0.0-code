"""Check completed campaign provenance, solver invariants and reported noise checks."""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'data/raw/controlled_20261117'
OUT = ROOT/'docs/research_refocus_20261117/campaign_audit.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    cells, seeds = [], []
    expected = {'baseline': {64: 16, 128: 16, 256: 16, 512: 16, 1024: 64},
                'period4': {128: 32, 256: 32}, 'period8': {128: 32, 256: 32},
                'modes_plus4': {128: 32, 256: 32}}
    for group, sizes in expected.items():
        for size, runs in sizes.items():
            path = RAW/group/f'theta0.25_L{size}.npz'
            with np.load(path, allow_pickle=False) as archive:
                z = {k: archive[k] for k in archive.files}
            spec = json.loads(str(z['spec']))
            assert spec['runs'] == runs and spec['L'] == size
            for name, file in [('source_sha256', 'controlled_roughening.py'),
                               ('solver_sha256', 'bd_core.py'), ('archive_sha256', 'control_archives.py')]:
                assert spec[name] == digest(ROOT/'code'/file)
            assert z['late_w2'].shape == (runs, 2, len(spec['observation_times']))
            assert spec['observation_times'][-1] <= spec['record_length']
            for k, array in z.items():
                if np.issubdtype(array.dtype, np.number):
                    assert np.all(np.isfinite(array)), (path, k)
            assert np.max(abs(z['profiles']-np.roll(z['profiles'], 1, axis=-1))) <= 1
            variance = z['profiles'].var(axis=-1)
            np.testing.assert_allclose(z['structure_factor'].sum(axis=-1)/size, variance, atol=1e-9)
            times = z['w_times'][0, 0]
            assert np.all(z['w_times'] == times)
            for j, stop in enumerate(spec['observation_times']):
                window = (times > .8*stop) & (times <= stop)
                np.testing.assert_array_equal(z['late_w2'][:, :, j], z['w_series'][:, :, window].mean(axis=-1))
            seeds.extend(map(int, z['variant_seeds'].ravel()))
            records = path.parent/(path.stem+'_runs')
            assert len(list(records.glob('run*.npz'))) == runs
            for i in range(runs):
                with np.load(records/f'run{i:06d}.npz', allow_pickle=False) as one:
                    for key in ['late_w2', 'profiles', 'w_series', 'variant_seeds']:
                        np.testing.assert_array_equal(one[key], z[key][i])
            cells.append(dict(file=path.relative_to(ROOT).as_posix(), sha256=digest(path),
                              group=group, L=size, runs=runs,
                              maximum_raw_correlation_error=np.max(abs(z['raw_correlation'].mean(axis=0)-z['theoretical_raw_correlation']), axis=-1).tolist(),
                              raw_variance_ratio=(z['raw_second_moment'].mean(axis=0)/z['theoretical_raw_covariance'][0]).tolist(),
                              maximum_binary_G_U_difference=float(np.max(abs(np.diff(z['sign_covariance'].mean(axis=0), axis=0)))),
                              deposition_mean=z['deposition_mean'].mean(axis=0).tolist(),
                              spatial_sign_offsets=z['spatial_sign_offsets'][0, 0].tolist(),
                              spatial_sign_correlations=z['spatial_sign_correlations'].mean(axis=0).tolist(),
                              measured_excess_kurtosis=z['raw_excess_kurtosis'].mean(axis=0).tolist(),
                              predicted_excess_kurtosis=z['predicted_excess_kurtosis'].tolist()))
    assert len(seeds) == len(set(seeds)) == 640
    calibration = RAW/'calibration_sensitivity'
    assert len(list(calibration.glob('fit_*.npz'))) == 65
    for path in calibration.glob('fit_*.npz'):
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z['metadata']))
            assert meta['producer_sha256'] == digest(ROOT/'code/calibration_sensitivity.py')
            assert meta['source_sha256'] == digest(ROOT/'data/raw/innovation_20261031_v3/theta0.25_L128.npz')
            assert meta['binary_sha256'] == digest(ROOT/'data/raw/innovation_binary_refined_20261031/theta0.25_L128.npz')
            assert z['late_w2'].shape == (16, 2) and np.all(z['late_w2'] > 0)
            assert np.all(np.isfinite(z['predicted_sign']))
    result = dict(passed=True, cells=cells, independent_interface_trajectories=len(seeds),
                  calibration='65 spectra including original; 16 common innovation records, not 1040 independent records.',
                  checks=['source/solver/archive hashes', 'complete run counts', 'unique variant seeds',
                          'finite values', 'single-step constraint', 'Parseval', 'late-window reaggregation',
                          'aggregate equals saved independent runs'],
                  noise_diagnostics='Descriptive deviations; no new post-hoc equivalence threshold imposed.')
    OUT.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('Campaign audit passed: 11 cells, 640 independent interface trajectories.')


if __name__ == '__main__':
    main()
