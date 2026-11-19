"""Whole-run inference for the independently specified time/cutoff pilot."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def paired_duration_ratio(w2, draws=2000, seed=20261117):
    rng = np.random.default_rng(seed)
    w2 = np.asarray(w2)
    boot = w2[rng.integers(len(w2), size=(draws, len(w2)))].mean(axis=1)
    means = w2.mean(axis=0)
    return dict(ratio=np.sqrt(means/means[0]).tolist(),
                ci95=np.quantile(np.sqrt(boot/boot[:, :1]), [.025, .975], axis=0).tolist())


def adjacent_contrast(sizes, gaussian, uniform, draws=2000, seed=20261117):
    rng = np.random.default_rng(seed)
    means, samples = [], []
    for records in [gaussian, uniform]:
        means.append(np.array([.5*np.log(w.mean(axis=0)) for w in records]))
        samples.append(np.array([.5*np.log(w[rng.integers(len(w), size=(draws, len(w)))].mean(axis=1))
                                 for w in records]))
    dx = np.diff(np.log(sizes))[:, None]
    alpha = [np.diff(m, axis=0)/dx for m in means]
    boot = [np.diff(b, axis=0)/dx[:, None] for b in samples]
    return dict(sizes=list(map(int, sizes)), alpha_gaussian=alpha[0].tolist(),
                alpha_uniform=alpha[1].tolist(), delta=(alpha[1]-alpha[0]).tolist(),
                alpha_ci95=[np.quantile(b, [.025, .975], axis=1).tolist() for b in boot],
                ci95=np.quantile(boot[1]-boot[0], [.025, .975], axis=1).tolist())


def load(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def analyse(root, out):
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261117)
    groups, sources = {}, []
    for path in sorted(root.glob('*/*.npz')):
        if path.parent.name not in ['baseline', 'period4', 'period8', 'modes_plus4']:
            continue
        if path.name.startswith('theta'):
            z = load(path)
            spec = json.loads(str(z['spec']))
            group = path.parent.name
            rows = groups.setdefault(group, [])
            w = z['late_w2']
            n = len(w)
            ix = rng.integers(n, size=(2000, n))
            row = dict(L=spec['L'], theta=spec['theta'], runs=n, T=spec['T'],
                       P=spec['record_length'], N=spec['N'], observation_times=spec['observation_times'],
                       width=np.sqrt(w.mean(axis=0)).tolist(),
                       width_ci95=np.quantile(np.sqrt(w[ix].mean(axis=1)), [.025, .975], axis=0).tolist(),
                       duration=[paired_duration_ratio(w[:, v]) for v in range(2)],
                       file=path.relative_to(ROOT).as_posix(),
                       lags=z['lags'].tolist(), four_spacings=z['four_spacings'].tolist())
            # Assess drift within the last observed window, not just T -> 2T.
            times = z['w_times'][0, 0]; stop = spec['observation_times'][-1]
            first = (times > .6*stop) & (times <= .8*stop)
            last = (times > .8*stop) & (times <= stop)
            row['tail_drift'] = [paired_duration_ratio(np.column_stack([
                z['w_series'][:, v, first].mean(axis=1),
                z['w_series'][:, v, last].mean(axis=1)])) for v in range(2)]
            for key in ['raw_correlation', 'sign_covariance', 'four_moment', 'four_cumulant',
                        'deposition_mean', 'spatial_sign_correlation']:
                row[key] = z[key].mean(axis=0).tolist()
                row[key+'_ci95'] = np.quantile(z[key][ix].mean(axis=1), [.025, .975], axis=0).tolist()
            row['theoretical_raw_correlation'] = z['theoretical_raw_correlation'].tolist()
            row['gaussian_sign_target'] = z['gaussian_sign_target'].tolist()
            row['predicted_excess_kurtosis'] = z['predicted_excess_kurtosis'].tolist()
            rows.append(row)
            sources.append(dict(file=row['file'], sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    # Existing 64-run small-size cells and new L1024 use the same stationary
    # period, mode-count, threshold and late-window definition. Never blend
    # cold-start production or shorter-duration pilot data into this curve.
    extension = None
    large = root/'baseline/theta0.25_L1024.npz'
    if large.exists():
        sizes = [64, 128, 256, 512, 1024]
        gs, us = [], []
        for size in sizes:
            path = large if size == 1024 else ROOT/f'data/raw/innovation_20261031_v3/theta0.25_L{size}.npz'
            z = load(path); spec = json.loads(str(z['spec']))
            assert spec['record_length'] == 2*spec['T']
            w = z['late_w2']
            gs.append(w[:, 0]); us.append(w[:, 1 if size == 1024 else 3])
            sources.append(dict(file=path.relative_to(ROOT).as_posix(), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        extension = adjacent_contrast(sizes, gs, us)
        extension['runs_per_size'] = list(map(len, gs))
        extension['protocol'] = 'Historical stationary same-kernel G/U for L64-512; new L1024; T(L),P=2T(L),N(2T).'
    comparisons = []
    for a, b in [('baseline', 'period4'), ('period4', 'period8'), ('period4', 'modes_plus4')]:
        for ca in groups.get(a, []):
            for cb in groups.get(b, []):
                if ca['L'] != cb['L']:
                    continue
                za, zb = load(ROOT/ca['file']), load(ROOT/cb['file'])
                for time in sorted(set(ca['observation_times']) & set(cb['observation_times'])):
                    ia, ib = ca['observation_times'].index(time), cb['observation_times'].index(time)
                    wa, wb = za['late_w2'][:, :, ia], zb['late_w2'][:, :, ib]
                    ba = wa[rng.integers(len(wa), size=(2000, len(wa)))].mean(axis=1)
                    bb = wb[rng.integers(len(wb), size=(2000, len(wb)))].mean(axis=1)
                    comparisons.append(dict(first=a, second=b, L=ca['L'], time=time,
                                            ratio=np.sqrt(wb.mean(axis=0)/wa.mean(axis=0)).tolist(),
                                            ci95=np.quantile(np.sqrt(bb/ba), [.025, .975], axis=0).tolist()))
    result = dict(groups=groups, extension=extension, comparisons=comparisons, sources=sources,
                  bootstrap=2000, seed=20261117,
                  inference='Whole independent runs; paired durations, independent ensembles/period groups; pointwise exploratory 95% intervals. Finite windows only.')
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('Analysed groups:', {k: [(r['L'], r['runs']) for r in v] for k, v in groups.items()}, flush=True)
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, default=ROOT/'data/raw/controlled_20261117')
    ap.add_argument('--out', type=Path, default=ROOT/'data/analysis/controlled_20261117')
    args = ap.parse_args()
    analyse(args.source.resolve(), args.out.resolve())
