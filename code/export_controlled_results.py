"""Export width ratios for the supplemental time and spectral controls."""
import json
from pathlib import Path

import numpy as np

from controlled_analysis import paired_duration_ratio

ROOT = Path(__file__).resolve().parents[1]


def interval_text(value, lo, hi, digits=4):
    return f'{value:.{digits}f}\\;[{lo:.{digits}f},{hi:.{digits}f}]'


def main():
    s = json.loads((ROOT/'data/analysis/controlled_20261117/summary.json').read_text())
    ext, groups = s['extension'], s['groups']
    assert ext and ext['runs_per_size'] == [64]*5, 'Large-size extension incomplete'
    assert sorted(c['L'] for c in groups['baseline']) == [64, 128, 256, 512, 1024]
    for name in ['period4', 'period8', 'modes_plus4']:
        assert sorted(c['L'] for c in groups[name]) == [128, 256]
        assert all(c['runs'] == 32 for c in groups[name])
    value = ext['delta'][-1][1]
    lo, hi = np.array(ext['ci95'])[:, -1, 1]
    final_contrast = interval_text(value, lo, hi)
    rows = []
    for size in [128, 256]:
        row = next(c for c in groups['period8'] if c['L'] == size)
        with np.load(ROOT/row['file']) as z:
            for v, name in enumerate(['G', 'U']):
                ratio = paired_duration_ratio(z['late_w2'][:, v, 1:3])
                rows.append(dict(control='time', L=size, input=name, estimate=ratio['ratio'][1],
                                 lo=ratio['ci95'][0][1], hi=ratio['ci95'][1][1]))
    for first, second, label in [('period4', 'period8', 'period'), ('period4', 'modes_plus4', 'modes')]:
        for c in s['comparisons']:
            base = next(x for x in groups.get(first, []) if x['L'] == c['L'])
            if c['first'] == first and c['second'] == second and c['time'] == 2*base['T']:
                for v, name in enumerate(['G', 'U']):
                    rows.append(dict(control=label, L=c['L'], input=name, estimate=c['ratio'][v],
                                     lo=c['ci95'][0][v], hi=c['ci95'][1][v]))
    calibration = json.loads((ROOT/'data/raw/controlled_20261117/calibration_sensitivity/summary.json').read_text())
    pct = calibration['width_ratio_percentiles']
    lines = [r'\newcommand{\ControlledWidthRows}{']
    for r in rows:
        row = (
            f"{r['control']} & {r['L']} & {r['input']} & {r['estimate']:.4f}"
            f" & [{r['lo']:.4f}, {r['hi']:.4f}]"
        )
        lines.append(row + r' \\')
    lines.append('}')
    target = ROOT/'docs/strict_binary_20261117/supplement_source/tex/tables/controlled_numbers.tex'
    target.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    (ROOT/'docs/research_refocus_20261117/exported_control_numbers.json').write_text(json.dumps(dict(
        final_adjacent_contrast=final_contrast, rows=rows, calibration_percentiles=pct), indent=2)+'\n')
    print('Exported supplemental control ratios.')


if __name__ == '__main__':
    main()
