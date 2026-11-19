import numpy as np

from controlled_analysis import adjacent_contrast, paired_duration_ratio


def test_adjacent_contrast_recovers_known_powers_with_unequal_ensemble_sizes():
    sizes = np.array([64, 128, 256])
    gaussian = [np.full((5, 2), float(n)) for n in sizes]
    uniform = [np.full((7, 2), n**1.6) for n in sizes]
    result = adjacent_contrast(sizes, gaussian, uniform, draws=50)
    np.testing.assert_allclose(result['alpha_gaussian'], .5)
    np.testing.assert_allclose(result['alpha_uniform'], .8)
    np.testing.assert_allclose(result['alpha_ci95'][0], .5)
    np.testing.assert_allclose(result['alpha_ci95'][1], .8)
    np.testing.assert_allclose(result['delta'], .3)
    np.testing.assert_allclose(result['ci95'], .3)


def test_duration_ratio_keeps_run_pairing():
    base = np.array([1., 3., 11., 70.])
    result = paired_duration_ratio(np.column_stack([base, 4*base]), draws=50)
    np.testing.assert_allclose(result['ratio'], [1., 2.])
    np.testing.assert_allclose(result['ci95'], [[1., 2.], [1., 2.]])
