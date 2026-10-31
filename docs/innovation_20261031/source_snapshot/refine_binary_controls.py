"""One documented refinement of the two failed binary covariance controls.

The tolerance remains 0.005. Fit a denser lag grid with 128 fresh 64-site
records; validate once with independent new 64-run, 64-site records.
Original failed archives are preserved. No validation-dependent fit retries.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy import fft
from numba import set_num_threads

from bd_core import ffgn_weights
from binary_covariance_control import fit_binary_spectrum
from control_archives import atomic_savez
from innovation_controls import ROOT, noise_field, noise_diagnostics, grow, gaussian_sign_covariance
from innovation_analysis import covariance_gate, analyze

SOURCE = ROOT/'data/raw/innovation_20261031_v3'
ORIGINAL = ROOT/'data/raw/innovation_binary_20261031_v3'
REFINED = ROOT/'data/raw/innovation_binary_refined_20261031'
SEED = 20261103


def all_site_diagnostics(eta, lags):
    chunks = [noise_diagnostics(eta[:, site:site+16], lags)
              for site in range(0, eta.shape[1], 16)]
    return {key: np.mean([c[key] for c in chunks], axis=0) for key in chunks[0]}


def fresh_seed(size, domain, run):
    return int(np.random.SeedSequence([SEED, size, domain, run]).generate_state(1)[0])


def training_record(run, spec, weights, rates, lags):
    eta = noise_field(weights, rates, spec['record_length'], 64,
                      fresh_seed(spec['L'], 0, run), 'uniform_shared')
    sign = 2*eta.T.astype(float)-1
    power = (abs(fft.rfft(sign, axis=-1))**2).mean(axis=0)
    return (fft.irfft(power, n=spec['record_length'])/spec['record_length'])[lags]


def validation_record(run, spec, amp, lags, weights, rates):
    set_num_threads(1)
    period = spec['record_length']
    rng = np.random.default_rng(fresh_seed(spec['L'], 1, run))
    z = rng.standard_normal((64, period))
    eta = (fft.irfft(fft.rfft(z, axis=-1)*amp, n=period, axis=-1) > 0).T.astype(np.uint8)
    result = {**grow(np.ascontiguousarray(eta[:, :spec['L']]), spec['T']),
              **all_site_diagnostics(eta, lags)}
    reference = noise_field(weights, rates, period, 64,
                            fresh_seed(spec['L'], 2, run), 'uniform_shared')
    diag = all_site_diagnostics(reference, lags)
    result['validation_reference_covariance'] = diag['sign_covariance']
    result['validation_reference_persistence'] = diag['persistence']
    return result


def refine(size):
    source = SOURCE/f'theta0.40_L{size}.npz'
    original = ORIGINAL/source.name
    with np.load(source, allow_pickle=False) as z:
        spec = json.loads(str(z['spec']))
    lags = np.unique(np.r_[np.arange(65), spec['lags']])
    weights, rates = ffgn_weights(.4, spec['N'])
    with ThreadPoolExecutor(max_workers=4) as pool:
        training = np.array(list(pool.map(
            lambda run: training_record(run, spec, weights, rates, lags), range(128))))
    spectrum, fit = fit_binary_spectrum(training.mean(axis=0), lags, spec['record_length'])
    rows = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i, row in enumerate(pool.map(
                lambda run: validation_record(run, spec, np.sqrt(spectrum), lags, weights, rates), range(64))):
            rows.append(row)
            if (i+1) % 16 == 0:
                print(f'Refined L={size}: {i+1}/64 fresh validation/interface runs', flush=True)
    data = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    gate = covariance_gate(data['sign_covariance'], data['validation_reference_covariance'])
    with np.load(original, allow_pickle=False) as z:
        previous_gate = covariance_gate(z['sign_covariance'], z['validation_reference_covariance'])
    metadata = dict(
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_file=source.relative_to(ROOT).as_posix(), spec=spec, fit=fit,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seed_root=SEED, training_runs=128, training_sites=64, validation_runs=64,
        validation_sites=64, covariance_equivalence_margin=.005,
        refinement='One denser-lag, larger-calibration refinement after initial failure; fresh independent validation; no threshold change.',
        previous_archive=original.relative_to(ROOT).as_posix(),
        previous_archive_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
        previous_gate=previous_gate, refined_gate=gate)
    path = REFINED/source.name
    atomic_savez(path, metadata=json.dumps(metadata), validation_lags=lags,
                 training_sign_covariance=training,
                 expected_sign_covariance=gaussian_sign_covariance(spectrum, spec['record_length'])[lags],
                 fitted_spectrum=spectrum, **data)
    print(f'L={size}: original bound={previous_gate["simultaneous_upper95"]:.6f}; '
          f'new bound={gate["simultaneous_upper95"]:.6f}; passed={gate["passed"]}', flush=True)
    return metadata


def prepare_subset():
    REFINED.mkdir(exist_ok=True)
    for size in [32, 64, 128, 256]:
        for theta in [.1, .25, .4]:
            name = f'theta{theta:.2f}_L{size}.npz'
            if theta != .4 or size not in [32, 64]:
                shutil.copyfile(ORIGINAL/name, REFINED/name)
    reports = [refine(size) for size in [32, 64]]
    report = ROOT/'docs/innovation_20261031/refinement_report.json'
    report.write_text(json.dumps(reports, indent=2), encoding='utf-8')
    return reports


def finish_full_when_ready():
    while not (ROOT/'docs/innovation_20261031/results/completion.json').exists():
        time.sleep(15)
    for theta in [.1, .25, .4]:
        name = f'theta{theta:.2f}_L512.npz'
        shutil.copyfile(ORIGINAL/name, REFINED/name)
    out = ROOT/'docs/innovation_20261031/results_refined'
    analyze(SOURCE, REFINED, out)
    from update_innovation_manuscript import export
    export(out/'summary.json')


if __name__ == '__main__':
    prepare_subset()
