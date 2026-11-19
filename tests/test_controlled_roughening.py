"""Regression checks for independent record, cutoff and observation controls."""
import importlib
import json

import numpy as np
import pytest
from scipy import fft, signal


def load_controls():
    # Fail explicitly on the missing feature before the new driver exists.
    assert importlib.util.find_spec('controlled_roughening') is not None, 'controlled driver missing'
    return importlib.import_module('controlled_roughening')


def test_growth_times_profiles_tail_windows_and_parseval():
    controls = load_controls()
    eta = np.random.default_rng(23).integers(0, 2, (73, 9), dtype=np.uint8)
    got = controls.growth_observables(eta, [17, 61], stride=2)
    h = np.zeros(9, dtype=np.int64)
    profiles, widths = [], []
    for t, e in enumerate(eta[:61], 1):
        h = np.maximum.reduce([h+e, np.roll(h, 1), np.roll(h, -1)])
        widths.append(h.var())
        if t in [17, 61]:
            profiles.append(h.copy())
    np.testing.assert_array_equal(got['w_times'], np.arange(2, 62, 2))
    np.testing.assert_allclose(got['w_series'], np.array(widths)[1::2], atol=2e-7)
    np.testing.assert_array_equal(got['profiles'], profiles)
    for j, duration in enumerate([17, 61]):
        mask = (got['w_times'] > .8*duration) & (got['w_times'] <= duration)
        assert got['late_w2'][j] == pytest.approx(got['w_series'][mask].mean())
    np.testing.assert_allclose(got['structure_factor'].sum(axis=-1)/9,
                               np.var(profiles, axis=-1), atol=1e-12)


@pytest.mark.parametrize('times', [[0], [74], [30, 20], [20, 20]])
def test_growth_rejects_invalid_duration_and_never_replays(times):
    controls = load_controls()
    with pytest.raises(ValueError):
        controls.growth_observables(np.ones((73, 4), dtype=np.uint8), times)


def test_complex_filter_equals_direct_periodic_ar_kernel():
    controls = load_controls()
    period = 127
    weights, rates = np.array([.7, .2]), np.array([.2, .91])
    transfer = controls.shared_transfer(weights, rates, period)
    z = np.random.default_rng(3).normal(size=(2, period))
    direct = np.zeros_like(z)
    for weight, rate in zip(weights, rates):
        ar = signal.lfilter([1], [1, -rate], z, axis=-1)
        ar += ar[:, -1, None]*rate**np.arange(1, period+1)/(1-rate**period)
        direct += weight*ar
    direct -= direct.mean(axis=-1, keepdims=True)
    np.testing.assert_allclose(fft.irfft(fft.rfft(z)*transfer, n=period), direct, atol=1e-13)
    assert np.max(abs(transfer.imag)) > .1
    assert transfer[0] == 0


def test_noise_statistics_equal_direct_periodic_moments():
    controls = load_controls()
    raw = np.random.default_rng(11).normal(size=(5, 63))
    lags, spacings = np.array([0, 1, 2, 3, 6, 9]), np.array([1, 3])
    got = controls.noise_statistics(raw, lags, spacings)
    signs = 2*(raw > 0).astype(float)-1
    raw_cov = np.array([(raw*np.roll(raw, -int(lag), axis=-1)).mean() for lag in lags])
    sign_cov = np.array([(signs*np.roll(signs, -int(lag), axis=-1)).mean() for lag in range(10)])
    np.testing.assert_allclose(got['raw_covariance'], raw_cov)
    np.testing.assert_allclose(got['raw_correlation'], raw_cov/raw_cov[0])
    np.testing.assert_allclose(got['sign_covariance'], sign_cov[lags])
    for j, d in enumerate(spacings):
        four = np.prod([np.roll(signs, -int(k*d), axis=-1) for k in range(4)], axis=0).mean()
        assert got['four_moment'][j] == pytest.approx(four)
        assert got['four_cumulant'][j] == pytest.approx(
            four-sign_cov[d]**2-sign_cov[2*d]**2-sign_cov[d]*sign_cov[3*d])
    assert got['deposition_mean'] == pytest.approx((raw > 0).mean())
    assert got['spatial_sign_correlation'] == pytest.approx((signs[:-1]*signs[1:]).mean())


def test_noise_field_has_unit_innovations_same_covariance_and_chunk_invariance():
    controls = load_controls()
    transfer = controls.shared_transfer(np.array([.8, .2]), np.array([.2, .8]), 1024)
    theory = fft.irfft(abs(transfer)**2, n=1024)
    for variant in ['gaussian_shared', 'uniform_shared']:
        eta, raw = controls.noise_field(transfer, 1024, 64, 110, variant,
                                        diagnostic_sites=64, rows_per_chunk=3)
        eta2, raw2 = controls.noise_field(transfer, 1024, 64, 110, variant,
                                          diagnostic_sites=64, rows_per_chunk=7)
        np.testing.assert_array_equal(eta, eta2)
        np.testing.assert_array_equal(raw, raw2)
        np.testing.assert_array_equal(eta.T, raw > 0)
        stats = controls.noise_statistics(raw, [0, 1, 2], [1])
        np.testing.assert_allclose(stats['raw_covariance'], theory[:3], atol=.025)
        assert abs(stats['spatial_sign_correlation']) < .025


def test_run_seeds_are_stable_and_distinct():
    controls = load_controls()
    values = [controls.run_seeds(81, .25, 8, i) for i in range(3)]
    assert values == [controls.run_seeds(81, .25, 8, i) for i in range(3)]
    assert len(set(np.array(values).ravel().tolist())) == 6


def test_cutoff_and_observation_do_not_change_mode_count(tmp_path):
    controls = load_controls()
    paths = [controls.run_cell(.25, 8, 1, 24, 71, 1, tmp_path/str(period),
                               period=period, observation_times=[12, 24])
             for period in [48, 96]]
    specs = [json.loads(str(np.load(p)['spec'])) for p in paths]
    assert specs[0]['N'] == specs[1]['N']
    assert specs[0]['observation_times'] == specs[1]['observation_times'] == [12, 24]
    assert [s['record_length'] for s in specs] == [48, 96]
    explicit = controls.run_cell(.25, 8, 1, 24, 71, 1, tmp_path/'fixed', n_modes=7)
    assert json.loads(str(np.load(explicit)['spec']))['N'] == 7


def test_atomic_resume_extends_runs_and_rejects_incompatible_provenance(tmp_path):
    controls = load_controls()
    path = controls.run_cell(.25, 8, 1, 24, 9, 1, tmp_path)
    run_path = tmp_path/'theta0.25_L8_runs'/'run000000.npz'
    original, modified = run_path.read_bytes(), run_path.stat().st_mtime_ns
    controls.run_cell(.25, 8, 2, 24, 9, 2, tmp_path)
    assert run_path.read_bytes() == original
    assert run_path.stat().st_mtime_ns == modified
    with np.load(path, allow_pickle=False) as z:
        assert z['late_w2'].shape == (2, 2, 2)
        assert z['profiles'].shape == (2, 2, 2, 8)
        assert z['raw_covariance'].shape[:2] == (2, 2)
        assert not np.array_equal(z['w_series'][0, 0], z['w_series'][0, 1])
        assert json.loads(str(z['spec']))['runs'] == 2
    with pytest.raises(ValueError, match='incompatible'):
        controls.run_cell(.25, 8, 2, 24, 10, 1, tmp_path)
    with np.load(run_path, allow_pickle=False) as z:
        archived = {k: z[k] for k in z.files}
    spec = json.loads(str(archived['spec']))
    spec['source_sha256'] = 'obsolete'
    archived['spec'] = json.dumps(spec)
    np.savez_compressed(run_path, **archived)
    with pytest.raises(ValueError, match='incompatible'):
        controls.run_cell(.25, 8, 2, 24, 9, 1, tmp_path)
