"""Evidence figures from frozen archives; no growth/noise simulation is performed.

Run with: python code/refocus_figures.py
All intervals are pointwise 95% whole-run bootstrap intervals (2,000 draws).
Binary-surrogate inference is conditional on its fixed fitted spectrum.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

from bd_core import ffgn_weights
from innovation_controls import spectral_filters
from innovation_analysis import four_sign_moment

ROOT = Path(__file__).resolve().parents[1]
SIZES = [32, 64, 128, 256, 512]
BOOT = 2000
SEED = 20261118
COLORS = {'G': '#336A91', 'U': '#B66A36', 'Gb': '#69528A'}
WIDTH = 183 / 25.4
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Times New Roman'],
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'mathtext.fontset': 'stix', 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.linewidth': .65,
    'xtick.major.width': .65, 'ytick.major.width': .65,
    'legend.frameon': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.facecolor': 'white', 'lines.linewidth': 1.1,
})


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean_interval(samples, rng, root=False):
    """Bootstrap whole runs; count weights avoid large indexed 3D copies."""
    samples = np.asarray(samples, dtype=float)
    assert np.all(np.isfinite(samples))
    n = len(samples)
    weights = rng.multinomial(n, np.full(n, 1 / n), size=BOOT) / n
    draws = weights @ samples
    mean = samples.mean(axis=0)
    if root:
        assert np.all(samples >= 0)
        mean, draws = np.sqrt(mean), np.sqrt(draws)
    return mean, np.quantile(draws, [.025, .975], axis=0)


def error(ax, x, y, ci, color, label=None, marker='o', ls='-', **kwargs):
    y, ci = np.asarray(y), np.asarray(ci)
    assert np.all(ci[0] <= y) and np.all(y <= ci[1])
    return ax.errorbar(x, y, yerr=[y - ci[0], ci[1] - y],
                      color=color, label=label, marker=marker,
                      ms=3.6, capsize=2, lw=1, elinewidth=.7,
                      ls=ls, **kwargs)


def panel(ax, letter, title):
    ax.set_title(title, loc='left', pad=10)
    ax.text(-.17, 1.065, letter, transform=ax.transAxes,
            fontsize=11, weight='bold', va='bottom')
    ax.tick_params(direction='out', length=3, pad=3)


def size_axis(ax, sizes=SIZES):
    ax.set_xscale('log', base=2)
    ax.set_xticks(sizes)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel('System size, L')


class EvidenceFigures:
    def __init__(self, out, notes):
        self.out, self.notes = out, notes
        out.mkdir(parents=True, exist_ok=True)
        notes.mkdir(parents=True, exist_ok=True)
        self.sources, self.rows, self.checks = {}, {}, []
        self.rng = np.random.default_rng(SEED)
        self.summary = self.read_json('docs/innovation_20261031/results_refined/summary.json')
        self.audit = self.read_json('docs/research_refocus_20261117/audit_results.json')
        self.old = self.read_json('data/analysis/revision_20261026/summary.json')
        self.archives, self.binary = {}, {}
        expected = {s['file']: s['sha256'] for s in self.summary['sources']}
        for size in SIZES:
            name = f'theta0.25_L{size}.npz'
            for directory, target in [('innovation_20261031_v3', self.archives),
                                      ('innovation_binary_refined_20261031', self.binary)]:
                rel = f'data/raw/{directory}/{name}'
                target[size] = self.read_npz(rel, expected[rel])
        assert self.audit['frozen_summary_exact_match']

    def register(self, rel, expected=None):
        path = ROOT / rel
        sha = digest(path)
        if expected is not None:
            assert sha == expected, f'Archive hash changed: {rel}'
        self.sources[rel] = {'file': rel, 'sha256': sha, 'checked_frozen_hash': expected is not None}
        return path

    def read_json(self, rel):
        return json.loads(self.register(rel).read_text(encoding='utf-8'))

    def read_npz(self, rel, expected=None):
        with np.load(self.register(rel, expected), allow_pickle=False) as z:
            return {key: z[key] for key in z.files}

    def add(self, figure, panel_id, series, x, y, ci=None, theta=.25, size=None, n=None):
        rows = self.rows.setdefault(figure, [])
        for i, (xi, yi) in enumerate(zip(np.atleast_1d(x), np.atleast_1d(y))):
            rows.append(dict(panel=panel_id, series=series, theta=theta, L=size,
                             x=float(xi), estimate=float(yi),
                             ci_low=None if ci is None else float(np.asarray(ci)[0, i]),
                             ci_high=None if ci is None else float(np.asarray(ci)[1, i]), n=n))

    def save(self, fig, name):
        for ext in ['pdf', 'svg', 'png']:
            fig.savefig(self.out / f'{name}.{ext}', dpi=400)
        plt.close(fig)
        with (self.out / f'{name}_source.csv').open('w', newline='', encoding='utf-8') as stream:
            rows = self.rows[name]
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def noise(self):
        name = 'fig1_noise_statistics'
        fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.35))
        fig.subplots_adjust(left=.10, right=.975, bottom=.20, top=.92, wspace=.34, hspace=.55)
        a, b, c, d = axes.flat
        z = self.archives[512]
        spec = json.loads(str(z['spec']))
        lags = np.asarray(spec['lags'])[1:]
        weights, rates = ffgn_weights(.25, spec['N'])
        transfer, _, _ = spectral_filters(tuple(weights), tuple(rates), spec['record_length'])
        raw_cov = np.fft.irfft(abs(transfer)**2, n=spec['record_length']) / 12
        rho = raw_cov[lags] / raw_cov[0]
        a.plot(lags, rho, color='.35', lw=1.2, ls='--', label='Common theory')
        self.add(name, 'a', 'common_theoretical', lags, rho, size=512)
        fresh = self.read_npz('data/raw/controlled_20261117/baseline/theta0.25_L512.npz')
        assert np.array_equal(fresh['lags'][1:], lags)
        for v, key, marker in [(0, 'G', 'o'), (1, 'U', 's')]:
            mean, ci = mean_interval(fresh['raw_correlation'][:, v, 1:], self.rng)
            error(a, lags, mean, ci, COLORS[key], key, marker, ls='none')
            self.add(name, 'a', key + '_measured', lags, mean, ci, size=512, n=len(fresh['raw_correlation']))
        a.set(xscale='log', xlabel=r'Temporal lag, $\tau$', ylabel=r'$\rho(\tau)=C(\tau)/C(0)$')
        a.legend(loc='upper right')
        panel(a, 'a', 'Continuous covariance: same kernel')
        for v, key, marker in [(0, 'G', 'o'), (3, 'U', 's')]:
            mean, ci = mean_interval(z['sign_covariance'][:, v, 1:], self.rng)
            error(b, lags, mean, ci, COLORS[key], key, marker)
            self.add(name, 'b', key, lags, mean, ci, size=512, n=64)
        b.set(xscale='log', xlabel=r'Temporal lag, $\tau$', ylabel=r'$C_s(\tau)=\langle s_t s_{t+\tau}\rangle$')
        b.legend(loc='upper right')
        panel(b, 'b', 'Binary drive: same kernel')
        cell = next(r for r in self.summary['cells'] if r['theta'] == .25 and r['L'] == 512)
        gate = cell['binary_covariance_gate']
        blags = np.asarray(cell['binary_validation_lags'])[1:]
        residual = -np.asarray(gate['residual'])[1:]
        ci = -np.asarray(gate['pointwise_ci95'])[::-1, 1:]
        error(c, blags, residual, ci, COLORS['Gb'], marker='o')
        c.axhline(0, color='.5', lw=.6)
        for boundary in [-.005, .005]:
            c.axhline(boundary, color='.5', ls='--', lw=.8)
        c.set(xscale='log', ylim=(-.0063, .0063), xlabel=r'Temporal lag, $\tau$',
              ylabel=r'$C_s(U)-C_s(G_b)$')
        c.text(.02, .985, 'Dashed: ±0.005 matching tolerance', transform=c.transAxes, va='top', fontsize=8)
        self.add(name, 'c', 'U_minus_Gb', blags, residual, ci, size=512, n=64)
        panel(c, 'c', 'Approximate binary match')
        selected = [next(r for r in self.summary['cells'] if r['theta'] == .25 and r['L'] == size) for size in SIZES]
        m4 = -np.array([r['four_sign_moment_difference'] for r in selected])
        m4ci = -np.array([r['four_sign_moment_difference_ci95'] for r in selected]).T[::-1]
        k4, k4ci = [], []
        for size in SIZES:
            record = self.binary[size]
            lags4 = record['validation_lags']
            ix = [int(np.flatnonzero(lags4 == k)[0]) for k in (1,2,3)]
            controls = []
            for prefix in ['', 'validation_reference_']:
                c = record['sign_covariance' if not prefix else prefix+'covariance']
                m = four_sign_moment(record[prefix+'persistence'], c, lags4)
                pooled = c[:, ix].mean(axis=0)
                draws = self.rng.integers(len(m), size=(BOOT,len(m)))
                def connected(m, c):
                    return m-c[...,0]**2-c[...,1]**2-c[...,0]*c[...,2]
                controls.append((connected(m.mean(),pooled),
                                 connected(m[draws].mean(axis=1),c[:,ix][draws].mean(axis=1))))
            k4.append(controls[1][0]-controls[0][0])
            k4ci.append(np.quantile(controls[1][1]-controls[0][1],[.025,.975]))
        k4ci = np.array(k4ci).T
        error(d, np.array(SIZES)*.97, m4, m4ci, COLORS['Gb'], r'$M_4$', 'o', ls='none')
        error(d, np.array(SIZES)*1.03, k4, k4ci, COLORS['G'], r'$K_4$', 's', ls='none')
        d.axhline(0, color='.5', lw=.6)
        d.set_ylabel(r'$U-G_b$')
        d.legend(loc='lower right')
        size_axis(d)
        panel(d, 'd', 'Four-time input contrasts')
        self.add(name, 'd', 'M4_U_minus_Gb', SIZES, m4, m4ci, n=64)
        self.add(name, 'd', 'K4_U_minus_Gb', SIZES, k4, k4ci, n=64)
        fig.text(.10, .067, 'θ = 0.25; a–c: L = 512. G/U: Gaussian/uniform innovations; $G_b$: binary-matched Gaussian.', fontsize=8)
        fig.text(.10, .027, 'Bars: pointwise 95% run-bootstrap; a: 16 fresh runs; b–d: 64 runs. Surrogate spectrum fixed.', fontsize=8)
        self.save(fig, name)

    def noise_validation(self):
        """Show the temporal and spatial checks for the assembled binary field."""
        name = 'figS1_noise_validation'
        noise = self.read_json('data/analysis/revision_20261026/noise_validation.json')
        fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, 2.75))
        fig.subplots_adjust(left=.105, right=.975, bottom=.23, top=.84, wspace=.37)
        for theta, color, marker in [(.25, COLORS['G'], 'o'), (.40, COLORS['U'], 's')]:
            row = next(r for r in noise['series'] if r['generator'] == 'fgn'
                       and r['theta'] == theta and r['L'] == 64)
            assert all(row['checks'][key] for key in
                       ('time_matches_finite_period_target', 'space_matches_zero'))
            for ax, letter, key, lags, target in [
                (a, 'a', 'time', noise['lags'], row['target_time']),
                (b, 'b', 'space', noise['spatial_lags'], np.zeros(len(noise['spatial_lags']))),
            ]:
                y = np.asarray(row[f'C_{key}'])
                se = np.asarray(row[f'C_{key}_err'])
                ax.errorbar(lags, y, yerr=se, color=color, marker=marker,
                            ms=3.6, capsize=2, elinewidth=.7, ls='none',
                            label=rf'$\theta={theta:.2f}$')
                self.add(name, letter, 'measured', lags, y, [y-se, y+se],
                         theta=theta, size=64, n=row['n_runs'])
                self.add(name, letter, 'prediction', lags, target, theta=theta, size=64)
            a.plot(noise['lags'], row['target_time'], ls='--', color=color, lw=.8, zorder=1)
        a.set_xscale('log', base=2)
        exponents = [0, 2, 4, 6]
        a.set_xticks([2**k for k in exponents], [rf'$2^{k}$' for k in exponents])
        a.set(xlabel=r'Temporal lag, $\Delta t$ (sweeps)', ylabel=r'$C_s(\Delta t,0)$')
        a.legend(loc='upper right')
        b.axhline(0, color='.45', ls='--', lw=.75)
        b.set(xticks=noise['spatial_lags'], xlabel=r'Spatial lag, $\Delta x$ (sites)',
              ylabel=r'$C_s(0,\Delta x)$')
        b.ticklabel_format(axis='y', style='plain')
        panel(a, 'a', 'Temporal correlations')
        panel(b, 'b', 'Spatial correlations')
        for row in self.rows[name]:
            row['sem_low'] = row.pop('ci_low')
            row['sem_high'] = row.pop('ci_high')
        self.save(fig, name)

    def widths(self):
        name = 'fig2_widths'
        fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 3.10))
        fig.subplots_adjust(left=.08, right=.98, top=.84, bottom=.30, wspace=.37)
        expected = {s['raw_file']: s['sha256'] for s in self.old['sources']}
        for ax, theta, letter in zip(axes[:2], [.25, .4], ['a', 'b']):
            for generator, label, color, marker in [
                ('shared', 'Shared AR', COLORS['U'], 's'),
                ('powerlaw', 'Power-law spectrum', COLORS['G'], 'o'),
                ('fgn', 'fGn spectrum', '#557D65', '^')]:
                means, intervals = [], []
                for size in SIZES:
                    row = next(r for r in self.old['cells'] if r['generator'] == generator and r['theta'] == theta and r['L'] == size)
                    z = self.read_npz(row['raw_file'], expected[row['raw_file']])
                    series = z['w_series'].astype(float)
                    count = max(1, int(series.shape[1] * .2))
                    samples = series[:, -count:].mean(axis=1)
                    mean, ci = mean_interval(samples, self.rng, root=True)
                    assert np.isclose(mean, row['W'], rtol=1e-13)
                    means.append(mean); intervals.append(ci)
                ci = np.array(intervals).T
                error(ax, SIZES, means, ci, color, label, marker)
                self.add(name, letter, generator, SIZES, means, ci, theta=theta, n=64)
            panel(ax, letter, f'Original protocol, θ = {theta:.2f}')
        large = self.read_npz('data/raw/controlled_20261117/baseline/theta0.25_L1024.npz')
        assert len(large['late_w2']) == 64
        for v, key, marker in [(0, 'G', 'o'), (3, 'U', 's')]:
            selected = [next(r for r in self.summary['cells'] if r['theta'] == .25 and r['L'] == size) for size in SIZES]
            means = np.array([r['widths'][v][1] for r in selected])
            ci = np.array([np.array(r['width_ci95'])[:, v, 1] for r in selected]).T
            last, last_ci = mean_interval(large['late_w2'][:, 0 if v == 0 else 1, 1], self.rng, root=True)
            means = np.r_[means, last]
            ci = np.column_stack([ci, last_ci])
            error(axes[2], SIZES+[1024], means, ci, COLORS[key], key, marker)
            self.add(name, 'c', key, SIZES+[1024], means, ci, n=64)
        panel(axes[2], 'c', 'Stationary input, θ = 0.25')
        for ax in axes:
            size_axis(ax)
            ax.set_yscale('log')
            ax.set_ylabel('Late-time width, W')
            ax.yaxis.set_major_formatter(ScalarFormatter())
        size_axis(axes[2], SIZES+[1024])
        axes[2].tick_params(axis='x', labelrotation=25)
        axes[0].legend(loc='upper left', fontsize=7.6, handlelength=1.4)
        axes[2].legend(loc='upper left')
        fig.text(.08, .105, 'a,b: original AR cold start / spectral period T. c: same-kernel G/U, period 2T; width at 2T.', fontsize=8)
        fig.text(.08, .045, r'Width = $\sqrt{\langle W^2\rangle}$ over final 20%. Bars: pointwise 95% bootstrap; n = 64 runs.', fontsize=8)
        self.save(fig, name)

    def flow(self):
        name = 'fig3_effective_flow'
        control = self.read_json('data/analysis/controlled_20261117/summary.json')
        extension = control['extension']
        assert extension and extension['runs_per_size'] == [64]*5
        for source in control['sources']:
            self.register(source['file'], source['sha256'])
        fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 3.8))
        fig.subplots_adjust(left=.08, right=.98, top=.84, bottom=.36, wspace=.40)
        positions = np.arange(4)
        ext_sizes = extension['sizes']
        for v, key, marker in [(0, 'G', 'o'), (1, 'U', 's')]:
            mean = np.array(extension['alpha_gaussian' if v == 0 else 'alpha_uniform'])[:, 1]
            ci = np.array(extension['alpha_ci95'])[v, :, :, 1]
            error(axes[0], positions + (v-.5)*.12, mean, ci, COLORS[key], key, marker, ls='none')
            self.add(name, 'a', key + '_2T', ext_sizes[:-1], mean, ci, n=64)
        axes[0].set_ylabel(r'$\alpha_{\mathrm{eff}}$')
        axes[0].legend(loc='lower left', ncol=2)
        panel(axes[0], 'a', 'Effective exponents at 2T')
        bounds = []
        for ax, key, letter, title in zip(axes[1:], ['same_kernel', 'binary_target'], ['b', 'c'],
                                          ['Same-kernel contrast', 'Approximate binary match']):
            if key == 'same_kernel':
                sizes = ext_sizes
                delta, intervals = np.array(extension['delta']), np.array(extension['ci95'])
            else:
                sizes = SIZES
                selected = [r for r in self.audit['local_indices'] if r['theta'] == .25 and r['comparison'] == key]
                assert [r['L'] for r in selected] == SIZES[:-1]
                delta = np.array([r['delta_U_minus_G'] for r in selected])
                intervals = np.stack([r['ci95'] for r in selected], axis=1)
            for duration, offset, marker, color in [(0, -.07, 'o', '#6D747B'), (1, .07, 's', COLORS['U'])]:
                means, ci = delta[:, duration], intervals[:, :, duration]
                bounds.extend(ci.ravel())
                error(ax, positions + offset, means, ci, color,
                      ['T(L)', '2T(L)'][duration], marker, ls='none')
                self.add(name, letter, key + '_' + ['T', '2T'][duration], sizes[:-1], means, ci,
                         n='G=64; U=64' if key == 'same_kernel' else 'Gb=64; U=32')
            ax.axhline(0, color='.5', lw=.8)
            ax.set_ylabel(r'$\Delta\alpha_{\mathrm{eff}}$')
            panel(ax, letter, title)
        limit = 1.15*max(.15, np.max(np.abs(bounds)))
        for ax in axes[1:]:
            ax.set_ylim(-limit, limit)
        for ax, sizes in zip(axes, [ext_sizes, ext_sizes, SIZES]):
            ax.set_xticks(positions, [f'{a}→{b}' for a, b in zip(sizes[:-1], sizes[1:])], rotation=25)
            ax.set_xlabel('Adjacent sizes, L → 2L')
        axes[1].legend(loc='upper left', ncol=2, fontsize=7.5, columnspacing=.5)
        fig.text(.08, .16,
                 r'$\theta=0.25$. (b) $\Delta\alpha_{\mathrm{eff}}=\alpha_{\mathrm{eff}}(U)-\alpha_{\mathrm{eff}}(G)$. '
                 r'(c) $\Delta\alpha_{\mathrm{eff}}=\alpha_{\mathrm{eff}}(U)-\alpha_{\mathrm{eff}}(G_b)$.', fontsize=8.5)
        fig.text(.08, .105, 'Bars: pointwise 95% run-bootstrap intervals. a,b: 64 runs/input; c: 32 held-out U and 64 $G_b$.', fontsize=8)
        fig.text(.08, .05, 'a,b extend to L = 1024; c ends at L = 512 and conditions on the fitted surrogate spectrum.', fontsize=8)
        self.save(fig, name)

    def structure(self):
        name = 'fig4_structure'
        fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 3.45), sharey=True)
        fig.subplots_adjust(left=.10, right=.98, top=.78, bottom=.29, wspace=.15)
        colors = ['#B8A68D', '#9A835F', '#708F8F', '#427083', '#253F64', '#252139']
        styles = [':', '--', '-.', '-', '-', '-']
        large = self.read_npz('data/raw/controlled_20261117/baseline/theta0.25_L1024.npz')
        max_parseval, max_fft = 0., 0.
        for ax, key, v, letter, title in zip(axes, ['G', 'U', 'Gb'], [0, 3, None], ['a', 'b', 'c'],
                                            ['Same-kernel G', 'Same-kernel U', 'Binary-matched $G_b$']):
            sizes = SIZES+[1024] if v is not None else SIZES
            for size, color, style in zip(sizes, colors, styles):
                if size == 1024:
                    column = 0 if key == 'G' else 1
                    standard = large['structure_factor'][:, column, 1]
                    profile = large['profiles'][:, column, 1].astype(float)
                else:
                    z = self.archives[size] if v is not None else self.binary[size]
                    stored = z['structure_factor'][:, v, 1] if v is not None else z['structure_factor'][:, 1]
                    profile = z['profiles'][:, v, 1].astype(float) if v is not None else z['profiles'][:, 1].astype(float)
                    standard = stored * size
                centered = profile - profile.mean(axis=1, keepdims=True)
                direct = abs(np.fft.fft(centered, axis=1))**2 / size
                assert np.allclose(standard, direct, rtol=1e-10, atol=1e-9)
                variance = (centered**2).mean(axis=1)
                relative = abs(standard.sum(axis=1) / size - variance) / variance
                assert np.max(relative) < 1e-12
                max_parseval = max(max_parseval, float(relative.max()))
                max_fft = max(max_fft, float(abs(standard - direct).max()))
                q = 2 * np.pi * np.arange(1, size // 2 + 1) / size
                mean, ci = mean_interval(standard[:, 1:size // 2 + 1], self.rng)
                assert np.all(q > 0) and np.all(ci > 0)
                ax.fill_between(q, ci[0], ci[1], color=color, alpha=.13, linewidth=0)
                ax.plot(q, mean, color=color, ls=style, lw=1, label=str(size))
                self.add(name, letter, key, q, mean, ci, size=size, n=64)
            ax.set(xscale='log', yscale='log', xlabel=r'Wave number, $q=2\pi m/L$')
            panel(ax, letter, title)
        axes[0].set_ylabel(r'$S(q)=\langle|\mathrm{FFT}(h-\bar h)|^2\rangle/L$')
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, ['L = ' + s for s in labels], loc='upper center', ncol=6,
                   bbox_to_anchor=(.54, .99), handlelength=2.6, columnspacing=1.4)
        fig.text(.10, .125, 'θ = 0.25; snapshot at 2T(L). Mean over 64 runs; bands: pointwise 95% run-bootstrap intervals.', fontsize=8)
        fig.text(.10, .065, 'All positive modes shown, including Nyquist. New L = 1024 data in a,b; binary-matched c ends at L = 512.', fontsize=8)
        self.checks.append(dict(check='structure_from_profiles_and_parseval', passed=True,
                               cells=17, runs_per_cell=64, max_relative_parseval_error=max_parseval,
                               max_absolute_fft_error=max_fft))
        self.save(fig, name)

    def manifest(self):
        result = dict(backend='Python/matplotlib', font='Times New Roman', width_mm=183,
                      bootstrap=BOOT, new_bootstrap_seed=SEED,
                      inference=('Pointwise run bootstrap for main comparisons with the fitted '
                                 'surrogate spectrum held fixed. Figure S1 uses one run standard error.'),
                      sources=list(self.sources.values()), checks=self.checks,
                      source_script={'file': 'code/refocus_figures.py', 'sha256': digest(Path(__file__))},
                      outputs=[{'file': str(p.resolve().relative_to(ROOT)).replace('\\', '/'), 'sha256': digest(p)}
                               for p in sorted(self.out.iterdir()) if p.suffix in {'.pdf', '.svg', '.png', '.csv'}])
        (self.notes / 'figure_provenance.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
        print(f'Exported {len(self.rows)} figures with {len(self.sources)} source files. Parseval checks passed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'figs/refocus_20261117')
    parser.add_argument('--notes', type=Path, default=ROOT / 'docs/research_refocus_20261117')
    args = parser.parse_args()
    figures = EvidenceFigures(args.out, args.notes)
    figures.noise()
    figures.widths()
    figures.flow()
    figures.structure()
    figures.noise_validation()
    figures.manifest()


if __name__ == '__main__':
    main()
