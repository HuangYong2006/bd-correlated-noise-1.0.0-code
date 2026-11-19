"""Two-innovation stationary controls with separate record, mode and time limits.

Unit-variance Gaussian and uniform innovations use the same complex circular AR
transfer with DC removed. Sites and variants are independent; durations share
prefixes of one record. Circular noise statistics never imply replay in growth.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import fft
from numba import set_num_threads

from bd_core import _simulate_from_field, ffgn_n_components, ffgn_weights, smooth_ceil
from control_archives import atomic_savez

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ('gaussian_shared', 'uniform_shared')
VERSION = 'controlled-circular-unitvariance-v1'


def shared_transfer(weights, rates, period):
    """Complex H=sum_n W_n/(1-r_n exp(-i omega)); only DC is removed."""
    phase = np.exp(-2j*np.pi*np.fft.rfftfreq(period))
    transfer = np.zeros_like(phase)
    for weight, rate in zip(weights, rates):
        transfer += weight/(1-rate*phase)
    transfer[0] = 0
    return transfer


def run_seeds(seed, theta, sites, run):
    sequence = np.random.SeedSequence([seed, int(round(theta*1e9)), sites, run, 20261117])
    return tuple(int(child.generate_state(1, dtype=np.uint64)[0])
                 for child in sequence.spawn(2))


def noise_field(transfer, period, sites, seed, variant, *, diagnostic_sites=16,
                rows_per_chunk=4):
    """Return binary (time,site) and the first diagnostic_sites continuous rows.

    Chunking changes neither random streams nor the time/site mapping. Only the
    binary full field is retained; continuous diagnostics are capped by caller.
    """
    if variant not in VARIANTS:
        raise ValueError(f'unknown variant: {variant}')
    rng = np.random.default_rng(seed)
    eta = np.empty((period, sites), dtype=np.uint8)
    raw = np.empty((min(sites, diagnostic_sites), period), dtype=np.float64)
    for start in range(0, sites, rows_per_chunk):
        n = min(rows_per_chunk, sites-start)
        z = (rng.standard_normal((n, period)) if variant == 'gaussian_shared'
             else rng.uniform(-np.sqrt(3), np.sqrt(3), (n, period)))
        spectrum = fft.rfft(z, axis=-1, workers=1)
        del z
        spectrum *= transfer
        x = fft.irfft(spectrum, n=period, axis=-1, workers=1)
        del spectrum
        eta[:, start:start+n] = (x > 0).T
        keep = min(n, max(0, len(raw)-start))
        if keep:
            raw[start:start+keep] = x[:keep]
        del x
    return eta, raw


def noise_statistics(raw, lags, four_spacings):
    """Pool known-zero-mean circular moments over sites (no record demeaning).

    The connected four-point estimator subtracts products of pooled population-
    mean-zero two-point estimators. It is a finite-sample plug-in estimator.
    Spatial checks use distinct physical sites without joining the end sites.
    """
    raw = np.asarray(raw, dtype=np.float64)
    sites, period = raw.shape
    lags = np.asarray(lags, dtype=int)
    four_spacings = np.asarray(four_spacings, dtype=int)
    if sites < 2 or np.any(lags < 0) or np.any(lags >= period):
        raise ValueError('need >=2 diagnostic sites and lags within record')
    if np.any(four_spacings <= 0) or np.any(3*four_spacings >= period):
        raise ValueError('four-point spacings must fit in record')
    required = np.unique(np.r_[0, lags, four_spacings, 2*four_spacings, 3*four_spacings])
    raw_cov = np.zeros(len(required))
    sign_cov = np.zeros(len(required))
    four = np.zeros(len(four_spacings))
    raw_sum = raw_second = raw_fourth = sign_sum = 0.
    spatial_offsets = np.array([d for d in [1, 2, 4] if d < sites])
    spatial = np.zeros(len(spatial_offsets))
    for i, row in enumerate(raw):
        spectrum = fft.rfft(row, workers=1)
        raw_cov += fft.irfft(abs(spectrum)**2, n=period, workers=1)[required]/period
        sign = 2*(row > 0).astype(float)-1
        spectrum = fft.rfft(sign, workers=1)
        sign_cov += fft.irfft(abs(spectrum)**2, n=period, workers=1)[required]/period
        raw_sum += row.mean()
        square = row*row
        raw_second += square.mean()
        raw_fourth += (square*square).mean()
        sign_sum += sign.mean()
        for j, d in enumerate(four_spacings):
            product = sign*np.roll(sign, -int(d))
            product *= np.roll(sign, -int(2*d))
            product *= np.roll(sign, -int(3*d))
            four[j] += product.mean()
        for j, d in enumerate(spatial_offsets):
            if i+d < sites:
                spatial[j] += (sign*(2*(raw[i+d] > 0)-1)).mean()/(sites-d)
    raw_cov /= sites
    sign_cov /= sites
    four /= sites
    c1 = sign_cov[np.searchsorted(required, four_spacings)]
    c2 = sign_cov[np.searchsorted(required, 2*four_spacings)]
    c3 = sign_cov[np.searchsorted(required, 3*four_spacings)]
    selected = np.searchsorted(required, lags)
    p = (sign_sum/sites+1)/2
    return dict(raw_covariance=raw_cov[selected], raw_correlation=raw_cov[selected]/raw_cov[0],
                sign_covariance=sign_cov[selected], four_moment=four,
                four_cumulant=four-c1*c1-c2*c2-c1*c3,
                raw_mean=raw_sum/sites, raw_second_moment=raw_second/sites,
                raw_excess_kurtosis=(raw_fourth/sites)/(raw_second/sites)**2-3,
                deposition_mean=p, diagnostic_deposition_mean=p,
                sign_mean=sign_sum/sites, spatial_sign_offsets=spatial_offsets,
                spatial_sign_correlations=spatial, spatial_sign_correlation=spatial[0])


def growth_observables(eta, observation_times, stride=1):
    """Evolve once to max(observation_times), retaining sampled W² and snapshots."""
    period, sites = eta.shape
    times = np.asarray(observation_times, dtype=np.int64)
    if (times.ndim != 1 or not len(times) or times[0] < 1 or times[-1] > period
            or np.any(np.diff(times) <= 0) or stride < 1):
        raise ValueError('observation times must increase within record; stride must be positive')
    w_times = np.arange(stride, int(times[-1])+1, stride, dtype=np.int64)
    windows = [(w_times > .8*t) & (w_times <= t) for t in times]
    if any(not mask.any() for mask in windows):
        raise ValueError('stride leaves an empty late-time window')
    profiles = np.zeros((1, len(times), sites), dtype=np.int32)
    widths = np.zeros((1, len(w_times)), dtype=np.float64)
    _simulate_from_field(eta[None], 0, sites, int(times[-1]), stride, times, profiles, widths)
    centered = profiles[0].astype(float)-profiles[0].mean(axis=-1, keepdims=True)
    return dict(w_times=w_times, w_series=widths[0], profiles=profiles[0],
                late_w2=np.array([widths[0, mask].mean() for mask in windows]),
                structure_factor=abs(fft.fft(centered, axis=-1, workers=1))**2/sites)


def _one_run(spec, run, transfer, directory):
    path = directory/f'run{run:06d}.npz'
    run_spec = {**spec, 'run': run}
    if path.exists():
        with np.load(path, allow_pickle=False) as archive:
            if json.loads(str(archive['spec'])) != run_spec:
                raise ValueError(f'incompatible existing run: {path}')
            return {k: archive[k] for k in archive.files if k != 'spec'}
    set_num_threads(1)
    start = time.perf_counter()
    seeds = run_seeds(spec['seed'], spec['theta'], spec['L'], run)
    rows = []
    for variant, seed in zip(VARIANTS, seeds):
        eta, raw = noise_field(transfer, spec['P'], spec['L'], seed, variant,
                               diagnostic_sites=spec['diagnostic_sites'],
                               rows_per_chunk=spec['rows_per_chunk'])
        diagnostics = noise_statistics(raw, spec['lags'], spec['four_spacings'])
        del raw
        diagnostics['deposition_mean'] = float(eta.mean())
        rows.append({**growth_observables(eta, spec['observation_times'], spec['stride']),
                     **diagnostics})
        del eta
    data = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    data['variant_seeds'] = np.array(seeds, dtype=np.uint64)
    data['wall_seconds'] = time.perf_counter()-start
    atomic_savez(path, spec=json.dumps(run_spec, sort_keys=True), **data)
    return data


def run_cell(theta, sites, runs, base_t, seed, workers, outdir, *, period=None,
             observation_times=None, n_modes=None, diagnostic_sites=16, rows_per_chunk=4):
    """Write resumable per-run archives and a (run,variant,...) aggregate.

    Default N depends on the baseline 2T only, never on the chosen P or t_obs.
    Increasing runs reuses completed independent runs without invalidating them.
    """
    period = 2*base_t if period is None else int(period)
    observation_times = [base_t, 2*base_t] if observation_times is None else list(observation_times)
    n_modes = ffgn_n_components(2*base_t) if n_modes is None else int(n_modes)
    if not (0 <= theta < .5) or min(sites, base_t, period, n_modes, runs) < 1:
        raise ValueError('invalid cell parameters')
    if workers not in [1, 2] or diagnostic_sites < 2 or sites < 2 or rows_per_chunk < 1:
        raise ValueError('use 1-2 workers, >=2 sites and a positive chunk size')
    if (not observation_times or observation_times[0] < 1 or observation_times[-1] > period
            or any(b <= a for a, b in zip(observation_times, observation_times[1:]))):
        raise ValueError('observation times must increase within record')
    lags = np.unique(np.r_[0, 1, 2, 4, 8, 16, 32, 64,
                           np.geomspace(128, max(128, base_t//2), 12).astype(int)])
    lags = lags[lags < period].tolist()
    spacings = [s for s in [1, 2, 4, 8, 16, 32, 64, 256] if 3*s < period]
    spec = dict(version=VERSION, theta=float(theta), L=sites, T=base_t, P=period,
                record_length=period, observation_times=observation_times, N=n_modes,
                seed=seed, variants=list(VARIANTS), lags=lags, four_spacings=spacings,
                diagnostic_sites=min(sites, diagnostic_sites), rows_per_chunk=rows_per_chunk,
                stride=max(1, min(base_t//2000, observation_times[0]//5)),
                innovation_variance=1., noise_covariance_estimator='circular-known-zero-mean',
                structure_factor_normalization='abs(FFT(h-mean(h)))**2/L',
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                solver_sha256=hashlib.sha256((ROOT/'code/bd_core.py').read_bytes()).hexdigest(),
                archive_sha256=hashlib.sha256((ROOT/'code/control_archives.py').read_bytes()).hexdigest())
    outdir = Path(outdir)
    stem = f'theta{theta:.2f}_L{sites}'
    directory = outdir/f'{stem}_runs'
    directory.mkdir(parents=True, exist_ok=True)
    weights, rates = ffgn_weights(theta, n_modes)
    transfer = shared_transfer(weights, rates, period)
    theory = fft.irfft(abs(transfer)**2, n=period, workers=1)
    kernel = fft.irfft(transfer, n=period, workers=1)
    uniform_excess = -1.2*np.sum(kernel**4)/np.sum(kernel**2)**2
    del kernel
    start = time.perf_counter()
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, result in enumerate(pool.map(lambda run: _one_run(spec, run, transfer, directory),
                                            range(runs)), 1):
            rows.append(result)
            if i % 8 == 0 or i == runs:
                print(f'theta={theta:.2f} L={sites} {i}/{runs} runs; '
                      f'{time.perf_counter()-start:.1f}s', flush=True)
    data = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    path = outdir/f'{stem}.npz'
    atomic_savez(path, spec=json.dumps({**spec, 'runs': runs}, sort_keys=True), **data,
                 lags=np.array(lags), four_spacings=np.array(spacings),
                 theoretical_raw_covariance=theory[lags],
                 theoretical_raw_correlation=theory[lags]/theory[0],
                 gaussian_sign_target=2/np.pi*np.arcsin(np.clip(theory[lags]/theory[0], -1, 1)),
                 predicted_excess_kurtosis=np.array([0., uniform_excess]), weights=weights, rates=rates)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--theta', type=float, default=.25)
    parser.add_argument('--sizes', type=int, nargs='+', default=[64, 128, 256, 512, 1024])
    parser.add_argument('--runs', type=int, default=16)
    parser.add_argument('--workers', type=int, choices=[1, 2], default=2)
    parser.add_argument('--seed', type=int, default=20261117)
    parser.add_argument('--out', type=Path, default=ROOT/'data/raw/controlled_20261117')
    parser.add_argument('--period-factor', type=float, default=2.)
    parser.add_argument('--observation-factors', type=float, nargs='+', default=[1., 2.])
    parser.add_argument('--n-modes', type=int)
    parser.add_argument('--base-t', type=int, help='override size-dependent baseline T')
    args = parser.parse_args()
    for sites in args.sizes:
        base_t = args.base_t if args.base_t is not None else smooth_ceil(round(15*sites**1.5))
        run_cell(args.theta, sites, args.runs, base_t, args.seed, args.workers, args.out,
                 period=round(args.period_factor*base_t),
                 observation_times=[round(f*base_t) for f in args.observation_factors],
                 n_modes=args.n_modes)


if __name__ == '__main__':
    main()
