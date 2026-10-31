"""Run-level inference and evidence report for the innovation campaign."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from innovation_controls import ROOT, VARIANTS

LABELS = ['G shared', 'G independent', 'G raw matched',
          'U shared', 'U independent', 'U raw matched', 'G binary matched']


def interval(values):
    return np.quantile(values, [.025, .975], axis=0).tolist()


def bootstrap_slopes(sizes, widths, n_boot=2000, seed=20261031):
    """One common run draw per size preserves coupling between variant columns."""
    rng = np.random.default_rng(seed)
    x = np.log(sizes); x -= x.mean()
    means, samples = [], []
    for w in widths:
        w = np.asarray(w)
        ix = rng.integers(len(w), size=(n_boot, len(w)))
        means.append(.5*np.log(w.mean(axis=0)))
        samples.append(.5*np.log(w[ix].mean(axis=1)))
    return np.tensordot(x, means, axes=1)/(x@x), np.tensordot(x, samples, axes=1)/(x@x)


def covariance_gate(surrogate, reference, margin=.005, n_boot=2000, seed=20261031):
    """95% simultaneous bootstrap upper bound on max absolute curve residual."""
    rng = np.random.default_rng(seed)
    surrogate, reference = np.asarray(surrogate), np.asarray(reference)
    a = surrogate[rng.integers(len(surrogate), size=(n_boot, len(surrogate)))].mean(axis=1)
    b = reference[rng.integers(len(reference), size=(n_boot, len(reference)))].mean(axis=1)
    residual = surrogate.mean(axis=0)-reference.mean(axis=0)
    upper = float(np.quantile(abs(a-b).max(axis=1), .95))
    return dict(max_absolute_residual=float(max(abs(residual))),
                simultaneous_upper95=upper, margin=margin, passed=upper < margin,
                residual=residual.tolist(), pointwise_ci95=interval(a-b))


def four_sign_moment(persistence, covariance, lags):
    """Consecutive four-sign moment using stationary sign-flip symmetry."""
    indices = np.searchsorted(lags, [1, 2, 3])
    return 16*persistence[:, 2]-1-covariance[:, indices] @ np.array([3., 2., 1.])


def covariance_validation_samples(data):
    if 'independent_validation_gaussian_covariance' in data:
        return (data['independent_validation_gaussian_covariance'],
                data['independent_validation_uniform_covariance'])
    return data['sign_covariance'], data['validation_reference_covariance']


def analyze(source, binary_source, out, partial=False, sizes=None):
    out.mkdir(parents=True, exist_ok=True)
    sizes = [32, 64, 128, 256, 512] if sizes is None else list(sizes)
    thetas = [.1, .25, .4]
    files = [source/f'theta{t:.2f}_L{size}.npz' for size in sizes for t in thetas]
    if not partial and any(not p.exists() for p in files):
        raise ValueError('formal campaign incomplete; use --partial only for monitoring')
    if not partial and any(not (binary_source/p.name).exists() for p in files):
        raise ValueError('binary campaign incomplete; use --partial only for monitoring')
    cells, provenance, loaded = [], [], {}
    rng = np.random.default_rng(20261031)
    for path in files:
        if not path.exists():
            continue
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        spec = json.loads(str(data['spec']))
        key = (spec['theta'], spec['L'])
        loaded[key] = data
        w = data['late_w2']
        ix = rng.integers(len(w), size=(2000, len(w)))
        mean_w = np.sqrt(w.mean(axis=0))
        duration = np.sqrt(w[:, :, 1].mean(axis=0)/w[:, :, 0].mean(axis=0))
        db = np.sqrt(w[ix, :, 1].mean(axis=1)/w[ix, :, 0].mean(axis=1))
        row = dict(theta=spec['theta'], L=spec['L'], runs=spec['runs'], T=spec['T'],
                   widths=mean_w.tolist(), duration_ratio=duration.tolist(), duration_ci95=interval(db),
                   width_ci95=interval(np.sqrt(w[ix].mean(axis=1))),
                   sign_covariance=data['sign_covariance'].mean(axis=0).tolist(), lags=spec['lags'],
                   deposition_mean=data['deposition_mean'].mean(axis=0).tolist(),
                   spatial_correlation=data['spatial_sign_correlation'].mean(axis=0).tolist(),
                   persistence=data['persistence'].mean(axis=0).tolist())
        uniform_diff = data['sign_covariance'][:, 3]-data['sign_covariance'][:, 5]
        row['uniform_raw_match_binary_max_residual'] = float(abs(uniform_diff.mean(axis=0)).max())
        gaussian_diff = data['sign_covariance'][:, 0]-data['sign_covariance'][:, 2]
        row['gaussian_positive_control_max_residual'] = float(abs(gaussian_diff.mean(axis=0)).max())
        binary = binary_source/path.name
        for p in [path, binary]:
            if p.exists():
                name = p.relative_to(ROOT) if p.is_relative_to(ROOT) else p.resolve()
                provenance.append(dict(file=name.as_posix(), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        if binary.exists():
            with np.load(binary, allow_pickle=False) as z:
                b = {k: z[k] for k in z.files}
            metadata = json.loads(str(b['metadata']))
            if metadata['source_sha256'] != hashlib.sha256(path.read_bytes()).hexdigest():
                raise ValueError(f'binary archive belongs to a different source: {binary}')
            if metadata['spec'] != spec:
                raise ValueError(f'binary archive parameters do not match source: {binary}')
            data['binary'] = b
            row['binary_covariance_gate'] = covariance_gate(*covariance_validation_samples(b))
            # Uniform held-out runs were not used to fit the surrogate spectrum.
            ref = w[32:, 3]
            bm = b['late_w2']
            ia = rng.integers(len(bm), size=(2000, len(bm)))
            ib = rng.integers(len(ref), size=(2000, len(ref)))
            ratio = np.sqrt(bm.mean(axis=0)/ref.mean(axis=0))
            rb = np.sqrt(bm[ia].mean(axis=1)/ref[ib].mean(axis=1))
            row['binary_matched_width_ratio'] = ratio.tolist()
            row['binary_matched_width_ratio_ci95'] = interval(rb)
            row['binary_widths'] = np.sqrt(bm.mean(axis=0)).tolist()
            row['binary_width_ci95'] = interval(np.sqrt(bm[ia].mean(axis=1)))
            row['binary_validation_lags'] = b['validation_lags'].tolist()
            row['binary_training_residual'] = metadata['fit']['max_training_residual']
            row['binary_calibration'] = {k: metadata[k] for k in
                ['training_runs', 'training_sites', 'validation_sites', 'refinement', 'previous_gate',
                 'previous_archive', 'additional_covariance_validation']
                if k in metadata}
            row['binary_persistence'] = b['persistence'].mean(axis=0).tolist()
            row['fresh_reference_persistence'] = b['validation_reference_persistence'].mean(axis=0).tolist()
            ip = rng.integers(len(bm), size=(2000, len(bm)))
            jp = rng.integers(len(bm), size=(2000, len(bm)))
            persistence = b['persistence']
            reference_persistence = b['validation_reference_persistence']
            row['persistence_difference'] = (persistence.mean(axis=0)-reference_persistence.mean(axis=0)).tolist()
            row['persistence_difference_ci95'] = interval(persistence[ip].mean(axis=1)-reference_persistence[jp].mean(axis=1))
            # Sign symmetry removes odd moments: 16 P(++++) = 1+3C1+2C2+C3+E[s0s1s2s3].
            fourth = four_sign_moment(persistence, b['sign_covariance'], b['validation_lags'])
            reference_fourth = four_sign_moment(reference_persistence, b['validation_reference_covariance'], b['validation_lags'])
            row['four_sign_moment_difference'] = float(fourth.mean()-reference_fourth.mean())
            row['four_sign_moment_difference_ci95'] = interval(fourth[ip].mean(axis=1)-reference_fourth[jp].mean(axis=1))
        cells.append(row)
    slopes = []
    comparisons = []
    pairs = [('gaussian_innovation', 0, 1), ('gaussian_raw_match', 0, 2),
             ('uniform_innovation', 3, 4), ('uniform_raw_match', 3, 5)]
    for theta in thetas:
        available = [s for s in sizes if (theta, s) in loaded]
        if len(available) < 3:
            continue
        inputs = [loaded[theta, s]['late_w2'] for s in available]
        central, boot = bootstrap_slopes(available, inputs)
        for v, name in enumerate(VARIANTS):
            slopes.append(dict(theta=theta, variant=name, sizes=available,
                               alpha=central[v].tolist(), ci95=interval(boot[:, v])))
        for name, a, b in pairs:
            delta, draws = central[b]-central[a], boot[:, b]-boot[:, a]
            comparisons.append(dict(theta=theta, comparison=name, sizes=available,
                                    delta_alpha=delta.tolist(), ci95=interval(draws)))
        if all('binary' in loaded[theta, s] for s in available):
            reference = [loaded[theta, s]['late_w2'][32:, 3] for s in available]
            surrogate = [loaded[theta, s]['binary']['late_w2'] for s in available]
            ac, ab = bootstrap_slopes(available, reference, seed=20261032)
            bc, bb = bootstrap_slopes(available, surrogate, seed=20261033)
            slopes.append(dict(theta=theta, variant='gaussian_binary_matched', sizes=available,
                               alpha=bc.tolist(), ci95=interval(bb)))
            relevant = [c for c in cells if c['theta'] == theta and c['L'] in available]
            comparisons.append(dict(theta=theta, comparison='binary_match_vs_heldout_uniform', sizes=available,
                                    delta_alpha=(bc-ac).tolist(), ci95=interval(bb-ab),
                                    all_covariance_gates_passed=all(c['binary_covariance_gate']['passed'] for c in relevant)))
    result = dict(status='PARTIAL' if partial else 'COMPLETE', requested_sizes=sizes,
                   cells=cells, slopes=slopes,
                   comparisons=comparisons, sources=provenance,
                   bootstrap=2000, inference='Whole runs; paired six-variant draws within size; independent sizes. Binary surrogate compared with held-out uniform runs; conditional on fixed fitted spectrum.',
                   limitations=['Finite record length=2T, DC removed; not the old cold-start protocol.',
                                'Finite-size slopes at T and 2T are not asymptotic exponents.',
                                'Binary covariance matching is approximate and checked only at sampled lags.',
                                'Bootstrap does not propagate spectrum-fit uncertainty or bound effects of unmeasured lag residuals.',
                                'Intervals are pointwise across physical comparisons; exploratory, not multiplicity-adjusted.'])
    (out/'summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    report(result, out)
    plot(result, loaded, out)
    print(f'{result["status"]}: {len(cells)}/{len(files)} cells; {sum("binary_covariance_gate" in c for c in cells)}/{len(files)} binary controls', flush=True)


def report(result, out):
    lines = ['# 噪声协方差与创新变量对照实验', '', f'状态：{result["status"]}。', '',
             '## 实验协议', '',
             '六类 AR 输入：高斯／均匀创新，各包含共享、独立、重配权重匹配原始协方差三种构造。高斯独立模态和匹配模态的标量谱因子采样在分布上严格等价于独立高斯模态之和。',
             '另拟合二值协方差匹配的高斯替代输入：基线方案使用前32次均匀共享运行的噪声统计作训练，采用非负 AR 谱混合。两组初始失败单元另作一次限定修正，以128次、每次64空间位置的新独立噪声记录加密滞后校准，再用64次、每次64位置的全新样本验证；阈值仍为0.005。若元数据包含additional_covariance_validation，则以保持拟合谱、界面数据不变的128次64位置独立纯噪声检验作为该单元的最终验证，原验证与高阶统计仍保留。界面比较保留后32次原均匀共享运行作为训练外参考。所有初始失败记录均保留。',
             '所有新输入均在长度2T的周期网格上去除DC，平坦界面从同一记录的前缀演化至T和2T。T沿用smooth_ceil(15 L^1.5)。这是一组控制初始化和截止的独立实验，不能直接当作旧冷启动数据的复现。',
             '宽度定义为sqrt(运行间平均的晚时段空间方差)。每个时段取该观测时长最后20%的规则采样点。跨尺寸拟合log(W)对log(L)得到有限尺寸有效斜率。', '',
             '## 有效斜率差：后者减前者', '',
             '| θ | 对照 | Δα(T) [95% CI] | Δα(2T) [95% CI] |',
             '|---|---|---|---|']
    for c in result['comparisons']:
        ci = np.asarray(c['ci95'])
        a = c['delta_alpha']
        lines.append(f'| {c["theta"]:.2f} | {c["comparison"]} | {a[0]:.5f} [{ci[0,0]:.5f}, {ci[1,0]:.5f}] | {a[1]:.5f} [{ci[0,1]:.5f}, {ci[1,1]:.5f}] |')
    lines += ['', '## 二值协方差的独立验证', '',
              '新生成的独立均匀共享参考和高斯替代样本，在0–64的每个整数滞后及原网格的长滞后上比较。下表为max|ΔCs|的同时bootstrap 95%上界；阈值为0.005。该阈值是本轮分析容差，不是普适的物理等效标准。未通过时，不能将界面残差解释为高阶统计效应。', '',
              '| θ | L | 最大残差 | 同时上界 | 通过 | W替代/W参考(2T) [95% CI] |',
              '|---|---|---|---|---|---|']
    for c in result['cells']:
        if 'binary_covariance_gate' not in c:
            continue
        g = c['binary_covariance_gate']; ci = np.array(c['binary_matched_width_ratio_ci95'])
        lines.append(f'| {c["theta"]:.2f} | {c["L"]} | {g["max_absolute_residual"]:.5f} | {g["simultaneous_upper95"]:.5f} | {g["passed"]} | {c["binary_matched_width_ratio"][1]:.4f} [{ci[0,1]:.4f}, {ci[1,1]:.4f}] |')
    lines += ['', '## 输入高阶统计检查', '',
              '利用整体正负翻转对称性和时间平稳性，由16P(++++)−1−3Cs(1)−2Cs(2)−Cs(3)估计连续四点符号矩。下表为高斯替代减独立均匀参考的差。区间是探索性的逐点运行bootstrap，不作多重比较校正；协方差检验未通过的单元不能用于隔离高阶效应。', '',
              '| θ | L | 四点矩差 [95% CI] |', '|---|---|---|']
    for c in result['cells']:
        if 'four_sign_moment_difference' in c:
            lo, hi = c['four_sign_moment_difference_ci95']
            lines.append(f'| {c["theta"]:.2f} | {c["L"]} | {c["four_sign_moment_difference"]:.5f} [{lo:.5f}, {hi:.5f}] |')
    lines += ['', '## 判读规则与结论边界', '',
              '共享／独立但权重相同的比较衡量创新耦合改变输入过程后的总效应。重配权重后，高斯过程的完整协方差和二值过程分布相同，应作为正对照；统计未显著不同不能自动证明等效。',
              '均匀输入的原始协方差匹配不保证二值相关匹配；两者界面差异不能直接归因于高阶时间统计。',
              '二值拟合通过独立验证时，可讨论在已测滞后和容差范围内，二阶相关匹配后仍存在的形态／宽度残差。但有限滞后匹配、谱拟合不确定度和剩余低频误差尚不能排除，故不能据此声称严格证明协方差不足、普适类变化或新临界点。',
              '2T/T晚时宽度比用于检查观测时长敏感性；偏离1说明该尺寸上的时长效应仍需保留。结构因子和高度差关联可辅助定位形态，持续沉积概率用于检查二值时间结构，均不单独构成机制证明。', '',
              '## 可复现性', '',
              '每个原始归档含逐运行数据、参数、种子和源代码哈希。summary.json列出输入文件SHA256。置信区间基于2000次独立运行bootstrap；二值替代与训练外参考独立抽样，条件于固定拟合谱，未重拟合每次bootstrap。',
              '完整运行和分析命令见 REPRODUCE.md。']
    (out/'report_zh.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def plot(result, loaded, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.spines.top': False,
                         'axes.spines.right': False, 'legend.frameon': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    colors = ['#366b99', '#d57a35', '#4b9b80', '#8c5b96', '#b64b4b', '#968337', '#292929']
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for j, theta in enumerate([.1, .25, .4]):
        available = sorted(s for th, s in loaded if th == theta)
        if not available:
            continue
        for v, label in enumerate(LABELS):
            values = []
            for size in available:
                z = loaded[theta, size]
                if v == 6 and 'binary' not in z:
                    break
                w = z['binary']['late_w2'][:, 1] if v == 6 else z['late_w2'][:, v, 1]
                values.append(np.sqrt(w.mean()))
            if len(values) == len(available):
                rows = [next(c for c in result['cells'] if c['theta']==theta and c['L']==size)
                        for size in available]
                bounds = [np.array(c['binary_width_ci95'])[:, 1] if v==6
                          else np.array(c['width_ci95'])[:, v, 1] for c in rows]
                bounds = np.array(bounds).T
                axes[0, j].errorbar(available, values,
                                     yerr=[np.array(values)-bounds[0], bounds[1]-np.array(values)],
                                     fmt='o-', color=colors[v], ms=3, lw=.8, label=label)
                axes[0, j].set(xscale='log', yscale='log')
        axes[0, j].set(title=f'θ={theta:.2f}, duration 2T', xlabel='L', ylabel='Late-time width')
        z = loaded[theta, max(available)]
        lag = np.array(json.loads(str(z['spec']))['lags'])
        for v in [0, 1, 3, 4, 5]:
            axes[1, j].plot(lag[1:], z['sign_covariance'][:, v].mean(axis=0)[1:],
                            'o-', ms=2, color=colors[v], label=LABELS[v])
        if 'binary' in z:
            b = z['binary']; lags = b['validation_lags']
            axes[1, j].plot(lags[1:], b['sign_covariance'].mean(axis=0)[1:], color=colors[6], label=LABELS[6])
        axes[1, j].set(xscale='log', xlabel='Temporal lag', ylabel='Sign covariance', title=f'Input check, L={max(available)}')
    axes[0, 0].legend(fontsize=6)
    axes[1, 0].legend(fontsize=6)
    fig.savefig(out/'innovation_controls.png', dpi=300)
    fig.savefig(out/'innovation_controls.pdf')
    fig.savefig(out/'innovation_controls.svg')
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    names = [('gaussian_innovation', 'G independent − shared'),
             ('gaussian_raw_match', 'G raw matched − shared'),
             ('uniform_innovation', 'U independent − shared'),
             ('uniform_raw_match', 'U raw matched − shared'),
             ('binary_match_vs_heldout_uniform', 'G binary matched − U held-out')]
    for j, theta in enumerate([.1, .25, .4]):
        for i, (key, label) in enumerate(names):
            rows = [c for c in result['comparisons'] if c['theta']==theta and c['comparison']==key]
            if not rows:
                continue
            c = rows[0]; ci = np.array(c['ci95'])
            for d, marker, offset in [(0, 'o', -.12), (1, 's', .12)]:
                x = c['delta_alpha'][d]
                axes[j].errorbar(x, i+offset, xerr=[[x-ci[0,d]], [ci[1,d]-x]],
                                 fmt=marker, ms=4, capsize=2, color=['#366b99', '#d57a35'][d],
                                 label=['T', '2T'][d] if i==0 else None)
        axes[j].axvline(0, color='.5', lw=.7)
        axes[j].set(yticks=range(5), yticklabels=[n[1] for n in names],
                    xlabel='Difference in finite-size slope (95% CI)', title=f'θ={theta:.2f}')
        axes[j].invert_yaxis()
    axes[0].legend()
    fig.savefig(out/'slope_differences.png', dpi=300)
    fig.savefig(out/'slope_differences.pdf')
    fig.savefig(out/'slope_differences.svg')
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
    for j, theta in enumerate([.1, .25, .4]):
        available = sorted(s for th, s in loaded if th == theta)
        if not available:
            continue
        size = max(available); z = loaded[theta, size]
        for v, label in enumerate(LABELS):
            if v == 6 and 'binary' not in z:
                continue
            sk = (z['binary']['structure_factor'][:, 1].mean(axis=0) if v == 6
                  else z['structure_factor'][:, v, 1].mean(axis=0))
            k = np.arange(1, size//2+1)
            axes[j].loglog(k, sk[k], '-', lw=1, color=colors[v], label=label)
        axes[j].set(title=f'θ={theta:.2f}, L={size}, 2T', xlabel='Mode k', ylabel='Structure factor')
    axes[0].legend(fontsize=6)
    fig.savefig(out/'structure_factors.png', dpi=300)
    fig.savefig(out/'structure_factors.pdf')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, default=ROOT/'data/raw/innovation_20261031_v3')
    ap.add_argument('--binary-source', type=Path, default=ROOT/'data/raw/innovation_binary_20261031_v3')
    ap.add_argument('--out', type=Path, default=ROOT/'docs/innovation_20261031/results')
    ap.add_argument('--partial', action='store_true')
    ap.add_argument('--sizes', nargs='+', type=int, default=[32, 64, 128, 256, 512])
    args = ap.parse_args()
    analyze(args.source.resolve(), args.binary_source.resolve(), args.out.resolve(), args.partial, args.sizes)


if __name__ == '__main__':
    main()
