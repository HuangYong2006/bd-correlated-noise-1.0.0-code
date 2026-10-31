"""One larger independent covariance-only check; fit and interfaces stay fixed."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import fft

from bd_core import ffgn_weights
from control_archives import atomic_savez
from innovation_analysis import covariance_gate
from innovation_controls import ROOT, noise_field


def sign_covariance(eta, lags):
    sign = 2*eta.T.astype(float)-1
    power = (abs(fft.rfft(sign, axis=-1))**2).mean(axis=0)
    return (fft.irfft(power, n=len(eta))/len(eta))[lags]


def validate():
    original = ROOT/'data/raw/innovation_binary_20261031_v3/theta0.40_L512.npz'
    target = ROOT/'data/raw/innovation_binary_refined_20261031/theta0.40_L512.npz'
    with np.load(original, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files}
    metadata = json.loads(str(data['metadata']))
    spec = metadata['spec']; period = spec['record_length']; lags = data['validation_lags']
    amp = np.sqrt(data['fitted_spectrum'])
    weights, rates = ffgn_weights(.4, spec['N'])

    def one(run):
        seeds = np.random.SeedSequence([20261105, 512, run]).spawn(2)
        rng = np.random.default_rng(seeds[0])
        eta = np.empty((period, 64), dtype=np.uint8)
        for start in range(0, 64, 8):
            white = rng.standard_normal((8, period))
            eta[:, start:start+8] = (fft.irfft(fft.rfft(white, axis=-1)*amp,
                                               n=period, axis=-1)>0).T
        gaussian = sign_covariance(eta, lags)
        reference = noise_field(weights, rates, period, 64,
                                int(seeds[1].generate_state(1)[0]), 'uniform_shared')
        return gaussian, sign_covariance(reference, lags)

    rows = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        for i, row in enumerate(pool.map(one, range(128))):
            rows.append(row)
            if (i+1) % 16 == 0:
                print(f'Fixed-spectrum L512 validation {i+1}/128', flush=True)
    gaussian, reference = np.asarray(rows).transpose(1, 0, 2)
    previous = covariance_gate(data['sign_covariance'], data['validation_reference_covariance'])
    gate = covariance_gate(gaussian, reference)
    metadata['additional_covariance_validation'] = dict(
        runs=128, sites=64, seed_root=20261105, fitted_spectrum_changed=False,
        interface_runs_changed=False, covariance_equivalence_margin=.005,
        previous_gate=previous, gate=gate,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        previous_archive=original.relative_to(ROOT).as_posix(),
        previous_archive_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
        protocol='One new independent covariance-only cohort after initial failure; no fit changes or further retries.')
    data['metadata'] = json.dumps(metadata)
    data['independent_validation_gaussian_covariance'] = gaussian
    data['independent_validation_uniform_covariance'] = reference
    atomic_savez(target, **data)
    report = ROOT/'docs/innovation_20261031/large_validation_report.json'
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(metadata['additional_covariance_validation'], indent=2), encoding='utf-8')
    print(f'Previous bound={previous["simultaneous_upper95"]:.6f}; new={gate["simultaneous_upper95"]:.6f}; passed={gate["passed"]}', flush=True)


if __name__ == '__main__':
    validate()
