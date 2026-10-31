import numpy as np
import pytest

import innovation_analysis as analysis


def test_slopes_average_squared_width_before_square_root():
    sizes = np.array([32, 64, 128, 256, 512])
    widths = [np.full((8, 2), size**1.4)*np.array([1., 4.]) for size in sizes]
    central, boot = analysis.bootstrap_slopes(sizes, widths, 100, 71)
    np.testing.assert_allclose(central, [.7, .7], atol=1e-12)
    np.testing.assert_allclose(boot, .7, atol=1e-12)


def test_covariance_gate_requires_entire_curve_within_margin():
    reference = np.ones((32, 5))*.2
    close = reference+.001
    far = reference.copy(); far[:, 3] += .02
    assert analysis.covariance_gate(close, reference, .005, 100)['passed']
    assert not analysis.covariance_gate(far, reference, .005, 100)['passed']


def test_extra_validation_does_not_replace_interface_or_moment_records():
    original = np.ones((64, 4))
    extra = np.zeros((128, 4))
    data = dict(sign_covariance=original, validation_reference_covariance=original,
                independent_validation_gaussian_covariance=extra,
                independent_validation_uniform_covariance=extra)
    a, b = analysis.covariance_validation_samples(data)
    assert a is extra and b is extra
    assert data['sign_covariance'] is original


def test_four_sign_moment_identity_on_symmetric_stationary_distribution():
    import itertools
    states = np.array(list(itertools.product([-1, 1], repeat=4)))
    weight = np.arange(1, 17, dtype=float)
    # Symmetrize under global sign flip and cyclic shifts.
    weight = (weight+weight[::-1])/2
    for _ in range(2):
        weight = sum(np.array([weight[np.flatnonzero((states == np.roll(s, shift)).all(axis=1))[0]]
                               for s in states]) for shift in range(4))/4
    weight /= weight.sum()
    covariance = np.array([[np.sum(weight*states[:, 0]*states[:, lag]) for lag in range(4)]])
    persistence = np.zeros((1, 6)); persistence[0, 2] = weight[(states == 1).all(axis=1)].sum()
    expected = np.sum(weight*states.prod(axis=1))
    np.testing.assert_allclose(analysis.four_sign_moment(persistence, covariance, np.arange(4)), expected, atol=1e-14)


def test_complete_report_requires_every_binary_control(tmp_path):
    source = tmp_path/'source'; source.mkdir()
    binary = tmp_path/'binary'; binary.mkdir()
    for size in [32, 64, 128, 256, 512]:
        for theta in [.1, .25, .4]:
            np.savez(source/f'theta{theta:.2f}_L{size}.npz', placeholder=0)
    with pytest.raises(ValueError, match='binary'):
        analysis.analyze(source, binary, tmp_path/'out')


def test_analysis_rejects_binary_control_from_different_raw_archive(tmp_path):
    import json
    source = tmp_path/'source'; source.mkdir()
    binary = tmp_path/'binary'; binary.mkdir()
    spec = dict(theta=.1, L=32, T=10, runs=2, lags=[0, 1])
    np.savez(source/'theta0.10_L32.npz', spec=json.dumps(spec),
             late_w2=np.ones((2, 6, 2)), sign_covariance=np.ones((2, 6, 2)),
             deposition_mean=np.ones((2, 6)), spatial_sign_correlation=np.zeros((2, 6)),
             persistence=np.zeros((2, 6, 6)))
    np.savez(binary/'theta0.10_L32.npz', metadata=json.dumps(dict(source_sha256='stale', spec=spec)))
    with pytest.raises(ValueError, match='different source'):
        analysis.analyze(source, binary, tmp_path/'out', partial=True)
