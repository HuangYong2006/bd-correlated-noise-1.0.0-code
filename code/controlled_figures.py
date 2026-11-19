"""Publication time/cutoff figure from completed independent control cohorts."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from controlled_analysis import paired_duration_ratio
from refocus_figures import COLORS, WIDTH, panel, error

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'figs/refocus_20261117'


def main():
    source = ROOT/'data/analysis/controlled_20261117/summary.json'
    summary = json.loads(source.read_text())
    groups = summary['groups']
    assert all(len(groups.get(g, [])) == 2 for g in ['period4', 'period8', 'modes_plus4'])
    rows, provenance = [], []
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.0))
    fig.subplots_adjust(left=.11, right=.98, top=.94, bottom=.19, hspace=.55, wspace=.35)
    a, b, c, d = axes.flat
    rng = np.random.default_rng(20261121)
    cell = next(x for x in groups['period8'] if x['L'] == 256)
    path = ROOT/cell['file']
    with np.load(path, allow_pickle=False) as archive:
        data = {k: archive[k] for k in archive.files}
    provenance.append({'file': cell['file'], 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    n = len(data['late_w2'])
    weights = rng.multinomial(n, np.full(n, 1/n), size=2000)/n
    # Every tenth stored point for display only; all samples enter late windows.
    select = np.arange(0, data['w_series'].shape[-1], 10)
    times = data['w_times'][0, 0, select]/cell['T']
    for v, key in enumerate(['G', 'U']):
        samples = data['w_series'][:, v, select]
        denominator = data['late_w2'][:, v, 0]
        mean = np.sqrt(samples.mean(axis=0)/denominator.mean())
        boot = np.sqrt((weights@samples)/(weights@denominator)[:, None])
        ci = np.quantile(boot, [.025, .975], axis=0)
        a.fill_between(times, ci[0], ci[1], color=COLORS[key], alpha=.15, lw=0)
        a.plot(times, mean, color=COLORS[key], label=key)
        for x, y, lo, hi in zip(times, mean, ci[0], ci[1]):
            rows.append(dict(panel='a', input=key, L=256, x=x, estimate=y, ci_low=lo, ci_high=hi))
    a.set(xscale='log', xlabel=r'Observation time, $t/T$', ylabel=r'$W(t)/\overline{W}(T)$')
    a.legend(loc='lower right', ncol=2)
    panel(a, 'a', 'Fixed record: P = 8T, L = 256')
    for v, key in enumerate(['G', 'U']):
        means, bounds = [], []
        for size in [128, 256]:
            row = next(x for x in groups['period8'] if x['L'] == size)
            path = ROOT/row['file']
            with np.load(path, allow_pickle=False) as archive:
                w = archive['late_w2'][:, v, 1:3]
            result = paired_duration_ratio(w)
            means.append(result['ratio'][1]); bounds.append(np.array(result['ci95'])[:, 1])
            rows.append(dict(panel='b', input=key, L=size, x=size, estimate=means[-1],
                             ci_low=bounds[-1][0], ci_high=bounds[-1][1]))
        offset = -.06 if v == 0 else .06
        error(b, np.arange(2)+offset, means, np.array(bounds).T, COLORS[key], key,
              marker=['o', 's'][v], ls='none')
    b.set(ylabel=r'$\overline{W}(4T)/\overline{W}(2T)$')
    panel(b, 'b', 'Duration sensitivity, P = 8T')
    for ax, first, second, letter, title, ylabel in [
        (c, 'period4', 'period8', 'c', 'Cutoff sensitivity, t = 2T', r'$\overline{W}(P=8T)/\overline{W}(P=4T)$'),
        (d, 'period4', 'modes_plus4', 'd', 'Mode truncation, P = 4T, t = 2T', r'$\overline{W}(N+4)/\overline{W}(N)$')]:
        for v, key in enumerate(['G', 'U']):
            means, bounds = [], []
            for size in [128, 256]:
                base = next(x for x in groups[first] if x['L'] == size)
                row = next(x for x in summary['comparisons'] if x['first'] == first and
                           x['second'] == second and x['L'] == size and x['time'] == 2*base['T'])
                means.append(row['ratio'][v]); bounds.append(np.array(row['ci95'])[:, v])
                rows.append(dict(panel=letter, input=key, L=size, x=size, estimate=means[-1],
                                 ci_low=bounds[-1][0], ci_high=bounds[-1][1]))
            error(ax, np.arange(2)+(-.06 if v == 0 else .06), means, np.array(bounds).T,
                  COLORS[key], marker=['o', 's'][v], ls='none')
        ax.set_ylabel(ylabel)
        panel(ax, letter, title)
    for ax in [b, c, d]:
        ax.axhline(1, color='.5', ls='--', lw=.8)
        ax.set(xticks=[0, 1], xticklabels=['128', '256'], xlabel='System size, L')
    fig.text(.11, .07, 'θ = 0.25; 32 runs per input/cell. Bars and bands: pointwise 95% whole-run bootstrap.', fontsize=8)
    fig.text(.11, .025, 'b: paired durations. c,d: independent cohorts. N is fixed when period or observation time changes.', fontsize=8)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT/f'fig5_convergence.{ext}', dpi=400)
    plt.close(fig)
    with (OUT/'fig5_convergence_source.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (OUT/'convergence_provenance.json').write_text(json.dumps(dict(
        source=str(source.relative_to(ROOT)), summary_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), sources=summary['sources'],
        plotted_time_sampling='Every tenth archived width sample; all samples retained in the archive and used in late-window estimates.'), indent=2)+'\n')
    print('Created fig5_convergence with source data and provenance.')


if __name__ == '__main__':
    main()
