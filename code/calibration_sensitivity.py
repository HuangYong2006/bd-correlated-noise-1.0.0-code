"""Exploratory fixed-L propagation of finite noise-calibration uncertainty.

Refit 64 run-bootstrap training samples. For each spectrum, use the same 16
fresh Gaussian innovation records to isolate its interface response. This is
a conditional sensitivity experiment, not an asymptotic equivalence test.
"""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import fft
from numba import set_num_threads

from binary_covariance_control import fit_binary_spectrum
from control_archives import atomic_savez
from innovation_controls import grow, gaussian_sign_covariance

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/raw/controlled_20261117/calibration_sensitivity'


def main():
    set_num_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    source = ROOT/'data/raw/innovation_20261031_v3/theta0.25_L128.npz'
    binary = ROOT/'data/raw/innovation_binary_refined_20261031/theta0.25_L128.npz'
    with np.load(source, allow_pickle=False) as z:
        spec = json.loads(str(z['spec']))
        training = z['sign_covariance'][:32, 3].copy()
    with np.load(binary, allow_pickle=False) as z:
        baseline = z['fitted_spectrum'].copy()
        validation_lags = z['validation_lags'].copy()
    period, sites = spec['record_length'], spec['L']
    rng = np.random.default_rng(2026111701)
    indices = rng.integers(32, size=(64, 32))
    producer_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    meta = dict(theta=.25, L=128, period=period, T=spec['T'],
                bootstrap_fits=64, interface_runs_per_fit=16, seed=2026111701,
                producer_sha256=producer_hash, source_sha256=source_hash,
                binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                interpretation='Conditional paired-interface response to bootstrap-refitted training spectra; '
                               'pointwise exploratory sensitivity, not a full physical equivalence bound.')
    rows = []
    start_time = time.perf_counter()
    for b in range(-1, 64):
        path = OUT/f'fit_{b+1:03d}.npz'
        if path.exists():
            with np.load(path, allow_pickle=False) as z:
                if json.loads(str(z['metadata'])) != meta:
                    raise ValueError(f'Incompatible cached result: {path}')
                rows.append({k: z[k] for k in z.files if k != 'metadata'})
            continue
        if b == -1:
            spectrum, coefficients = baseline, np.empty(0)
        else:
            spectrum, fit = fit_binary_spectrum(training[indices[b]].mean(axis=0), spec['lags'], period)
            coefficients = np.asarray(fit['coefficients'])
        assert spectrum[0] == 0 and np.all(spectrum >= 0)
        amp = np.sqrt(spectrum)
        widths = []
        for run in range(16):
            # Common random records across refits, independent of all old data.
            r = np.random.default_rng(np.random.SeedSequence([2026111702, run]))
            eta = np.empty((period, sites), dtype=np.uint8)
            for first in range(0, sites, 8):
                z = r.standard_normal((min(8, sites-first), period))
                eta[:, first:first+len(z)] = (fft.irfft(fft.rfft(z, axis=-1)*amp, n=period, axis=-1)>0).T
            widths.append(grow(eta, spec['T'])['late_w2'])
        row = dict(late_w2=np.asarray(widths), coefficients=coefficients,
                   predicted_sign=gaussian_sign_covariance(spectrum, period)[validation_lags])
        atomic_savez(path, metadata=json.dumps(meta), **row)
        rows.append(row)
        if (b+1) % 8 == 0:
            print(f'Calibration sensitivity {b+1}/64 fits, {time.perf_counter()-start_time:.1f}s', flush=True)
    w = np.array([r['late_w2'] for r in rows])
    sign = np.array([r['predicted_sign'] for r in rows])
    ratio = np.sqrt(w[1:].mean(axis=1)/w[0].mean(axis=0))
    residual = np.max(abs(sign[1:]-sign[0]), axis=1)
    summary = dict(metadata=meta, validation_lags=validation_lags.tolist(),
                   baseline_width=np.sqrt(w[0].mean(axis=0)).tolist(),
                   width_ratio_per_refit=ratio.tolist(),
                   width_ratio_percentiles=np.quantile(ratio, [.025, .5, .975], axis=0).tolist(),
                   width_ratio_range=[ratio.min(axis=0).tolist(), ratio.max(axis=0).tolist()],
                   max_sign_change_per_refit=residual.tolist(),
                   sign_change_percentiles=np.quantile(residual, [.025, .5, .975]).tolist())
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print('Completed conditional calibration sensitivity', flush=True)


if __name__ == '__main__':
    main()
