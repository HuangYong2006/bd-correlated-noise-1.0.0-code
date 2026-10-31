"""Independent NumPy/SciPy audit of an extracted supplement.

No project generator or analysis functions are imported. This checks archived
failures, current gates, regression/interval estimates, recalibration spectra
and one regenerated calibration record per refined cell.
"""
import hashlib
import json
from math import gamma
from pathlib import Path

import numpy as np


def load(path):
    with np.load(path, allow_pickle=False) as archive:
        return {k: archive[k] for k in archive.files}


def upper_bound(a, b):
    rng = np.random.default_rng(20261031)
    aa = a[rng.integers(len(a), size=(2000, len(a)))].mean(axis=1)
    bb = b[rng.integers(len(b), size=(2000, len(b)))].mean(axis=1)
    return float(np.quantile(np.max(np.abs(aa-bb), axis=1), .95))


def regress(sizes, records, seed):
    rng = np.random.default_rng(seed)
    y, boot = [], []
    for a in records:
        y.append(.5*np.log(a.mean(axis=0)))
        indices = rng.integers(len(a), size=(2000, len(a)))
        boot.append(.5*np.log(a[indices].mean(axis=1)))
    design = np.column_stack([np.log(sizes), np.ones(len(sizes))])
    y, boot = np.array(y), np.array(boot)
    slope = np.linalg.lstsq(design, y.reshape(len(sizes), -1), rcond=None)[0][0]
    draws = np.linalg.lstsq(design, boot.reshape(len(sizes), -1), rcond=None)[0][0]
    return slope.reshape(y.shape[1:]), draws.reshape(boot.shape[1:])


def audit(root):
    manifest = json.loads((root/'MANIFEST.json').read_text())
    for item in manifest['entries']:
        assert hashlib.sha256((root/item['file']).read_bytes()).hexdigest() == item['sha256']
    summary = json.loads((root/'summary.json').read_text())
    seeds = json.loads((root/'SEEDS.json').read_text())
    sizes = manifest['sizes']
    original_folder = root/'data/raw/innovation_binary_20261031_v3'
    chosen_folder = root/'data/raw/innovation_binary_refined_20261031'
    source_folder = root/'data/raw/innovation_20261031_v3'
    failures, checks = [], []
    for size in sizes:
        for theta in [.1, .25, .4]:
            name = f'theta{theta:.2f}_L{size}.npz'
            old, new = load(original_folder/name), load(chosen_folder/name)
            previous = upper_bound(old['sign_covariance'], old['validation_reference_covariance'])
            if previous >= .005:
                failures.append(dict(theta=theta, L=size, original_upper95=previous,
                                     file=(original_folder/name).relative_to(root).as_posix()))
            extra = 'independent_validation_gaussian_covariance' in new
            a = new['independent_validation_gaussian_covariance'] if extra else new['sign_covariance']
            b = new['independent_validation_uniform_covariance'] if extra else new['validation_reference_covariance']
            bound = upper_bound(a, b)
            cell = next(c for c in summary['cells'] if c['theta']==theta and c['L']==size)
            np.testing.assert_allclose(bound, cell['binary_covariance_gate']['simultaneous_upper95'], atol=1e-14)
            lag = new['validation_lags']
            i1, i2, i3 = [np.flatnonzero(lag==i)[0] for i in [1, 2, 3]]
            def fourth(p, c):
                return 16*p[:, 2]-1-3*c[:, i1]-2*c[:, i2]-c[:, i3]
            delta = fourth(new['persistence'], new['sign_covariance']).mean()
            delta -= fourth(new['validation_reference_persistence'],
                            new['validation_reference_covariance']).mean()
            np.testing.assert_allclose(delta, cell['four_sign_moment_difference'], atol=1e-12)
            if extra:
                for key in ['fitted_spectrum', 'late_w2', 'profiles', 'sign_covariance',
                            'validation_reference_covariance']:
                    np.testing.assert_array_equal(new[key], old[key])
            checks.append(dict(theta=theta, L=size, current_upper95=bound,
                               passed=bound < .005, independent_covariance_only=extra))
    for theta in [.1, .25, .4]:
        records = [load(source_folder/f'theta{theta:.2f}_L{s}.npz')['late_w2'] for s in sizes]
        central, draws = regress(sizes, records, 20261031)
        variants = ['gaussian_shared', 'gaussian_independent', 'gaussian_matched',
                    'uniform_shared', 'uniform_independent', 'uniform_matched']
        for i, variant in enumerate(variants):
            row = next(r for r in summary['slopes'] if r['theta']==theta and r['variant']==variant)
            np.testing.assert_allclose(central[i], row['alpha'], atol=1e-12)
            np.testing.assert_allclose(np.quantile(draws[:, i], [.025, .975], axis=0), row['ci95'], atol=1e-12)
        reference = [a[32:, 3] for a in records]
        target = [load(chosen_folder/f'theta{theta:.2f}_L{s}.npz')['late_w2'] for s in sizes]
        ac, ab = regress(sizes, reference, 20261032)
        bc, bb = regress(sizes, target, 20261033)
        row = next(r for r in summary['comparisons'] if r['theta']==theta
                   and r['comparison']=='binary_match_vs_heldout_uniform')
        np.testing.assert_allclose(bc-ac, row['delta_alpha'], atol=1e-12)
        np.testing.assert_allclose(np.quantile(bb-ab, [.025, .975], axis=0), row['ci95'], atol=1e-12)
    calibration_checks = []
    for size in [32, 64]:
        z = load(chosen_folder/f'theta0.40_L{size}.npz')
        meta = json.loads(str(z['metadata'])); spec = meta['spec']; period = spec['record_length']
        rates = np.r_[0., np.exp(-np.geomspace(6, 1/period, 60))]
        phase = np.exp(-2j*np.pi*np.fft.rfftfreq(period))
        basis = 1/abs(1-rates[:, None]*phase)**2
        basis[:, 0] = 0
        basis /= np.fft.irfft(basis, n=period, axis=-1)[:, 0, None]
        reconstructed = np.array(meta['fit']['coefficients']) @ basis
        np.testing.assert_allclose(reconstructed, z['fitted_spectrum'], rtol=1e-11, atol=1e-12)
        seed_row = next(r for r in seeds['recalibration'] if r['L']==size)
        seed = seed_row['training_uniform_seed32'][0]
        assert seed == int(np.random.SeedSequence([20261103, size, 0, 0]).generate_state(1)[0])
        rng = np.random.default_rng(seed)
        u = 6*2.**(-np.arange(1, spec['N']+1))
        r = np.exp(-u)
        k = (2**(.5-.4)-2**(.4-.5))/gamma(2-2*.4)
        w = np.sqrt((1-r*r)*k*u**(1-2*.4))
        transfer = sum(weight/(1-rate*phase) for weight, rate in zip(w, r))
        transfer[0] = 0
        sign = []
        for _ in range(8):
            white = rng.uniform(-.5, .5, (8, period))
            eta = np.fft.irfft(np.fft.rfft(white, axis=-1)*transfer, n=period, axis=-1)
            sign.append(2*(eta > 0).astype(float)-1)
        sign = np.concatenate(sign)
        power = (abs(np.fft.rfft(sign, axis=-1))**2).mean(axis=0)
        covariance = np.fft.irfft(power, n=period)/period
        err = float(np.max(abs(covariance[z['validation_lags']]-z['training_sign_covariance'][0])))
        assert err < 1e-12
        calibration_checks.append(dict(L=size, first_training_seed=seed,
                                       regenerated_training_max_error=err,
                                       fitted_spectrum_reconstructed=True))
    return dict(payload_hashes_verified=len(manifest['entries']), original_failures=failures,
                final_validation=checks, all_final_checks_pass=all(c['passed'] for c in checks),
                calibration_checks=calibration_checks,
                slope_and_bootstrap_intervals='independently reproduced by least-squares regression',
                four_sign_moments='all independently reproduced',
                imports_project_analysis_code=False)


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    report = audit(root)
    (root/'INDEPENDENT_AUDIT.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))
