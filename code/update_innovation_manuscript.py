"""Export a fixed control summary into the article, figure and compiled PDF."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from submission_archive import build_data_archive, build_submission_bundle

from revision_figures import configure_style, canvas, save
from innovation_controls import ROOT


def export(summary, figure_only=False):
    result = json.loads(summary.read_text(encoding='utf-8'))
    comparisons = result['comparisons']
    sizes = sorted(set.intersection(*[{c['L'] for c in result['cells'] if c['theta']==theta}
                                     for theta in [.1, .25, .4]]))
    if len(sizes) not in [4, 5] or any(c['sizes'] != sizes for c in comparisons):
        raise ValueError('need the same four or five completed sizes at every theta')
    figdir = ROOT/'figs/innovation_20261031'
    figdir.mkdir(parents=True, exist_ok=True)
    configure_style()
    fig, axes = canvas(height=99)
    fig.subplots_adjust(bottom=.18, top=.89, wspace=.42)
    for key, label, color, marker, offset in [
            ('gaussian_innovation', 'Gaussian: independent', '#1f6fb4', 'o', -.006),
            ('uniform_innovation', 'Uniform: independent', '#d1495b', 's', -.002),
            ('gaussian_raw_match', 'Gaussian: raw matched', '#4d7c69', '^', .002),
            ('uniform_raw_match', 'Uniform: raw matched', '#807565', 'D', .006)]:
        rows = [next(c for c in comparisons if c['theta']==theta and c['comparison']==key)
                for theta in [.1, .25, .4]]
        y = np.array([c['delta_alpha'][1] for c in rows])
        low = np.array([c['ci95'][0][1] for c in rows])
        high = np.array([c['ci95'][1][1] for c in rows])
        axes[0].errorbar(np.array([.1, .25, .4])+offset, y, yerr=[y-low, high-y],
                         fmt=marker+'-', color=color, ms=3, capsize=2, lw=.8, label=label)
    axes[0].axhline(0, color='.5', lw=.6)
    axes[0].set(xlabel='Temporal correlation parameter, θ',
                ylabel='Slope difference relative to shared input', xticks=[.1, .25, .4],
                ylim=(-.125, .035), xlim=(.07, .43))
    axes[0].legend(loc='lower left', fontsize=6.5)
    cells = [next(c for c in result['cells'] if c['theta']==theta and c['L']==256)
             for theta in [.1, .25, .4]]
    y = np.array([c['four_sign_moment_difference'] for c in cells])
    low = np.array([c['four_sign_moment_difference_ci95'][0] for c in cells])
    high = np.array([c['four_sign_moment_difference_ci95'][1] for c in cells])
    axes[1].errorbar([.1, .25, .4], y, yerr=[y-low, high-y],
                     fmt='o-', color='#443e68', ms=3.5, capsize=2)
    axes[1].axhline(0, color='.5', lw=.6)
    axes[1].set(xlabel='Temporal correlation parameter, θ',
                ylabel='Four-sign moment: Gaussian target − uniform',
                xticks=[.1, .25, .4], xlim=(.07, .43))
    outputs = save(fig, figdir, 'fig7_innovation_controls')
    if figure_only:
        print(f'Reproduced paper Figure 7: {figdir}', flush=True)
        return
    failed = [c for c in result['cells'] if not c['binary_covariance_gate']['passed']]
    failed_text = ', '.join(r'$\theta='+str(c['theta'])+r'$, $L='+str(c['L'])+'$' for c in failed)
    raw = [c['uniform_raw_match_binary_max_residual'] for c in result['cells']]
    binary_rows = [c for c in comparisons if c['comparison']=='binary_match_vs_heldout_uniform']
    eligible = {theta: all(c['binary_covariance_gate']['passed'] for c in result['cells']
                          if c['theta']==theta) for theta in [.1, .25, .4]}
    validated = [c for c in binary_rows if eligible[c['theta']]]
    unresolved = all(c['ci95'][0][1] <= 0 <= c['ci95'][1][1] for c in validated)
    stage = (r'The additional grid is restricted to the four completed sizes $L=32,64,128,256$; larger-size extensions are excluded from these estimates.'
             if len(sizes)==4 else r'The additional grid uses all five sizes $L=32,64,128,256,512$.')
    finding = ('Among parameter grids passing every covariance check, the surrogate slope contrasts at $2T$ are unresolved within the pointwise intervals in Table~\\ref{tab:innovation}.'
               if unresolved else 'Table~\\ref{tab:innovation} reports the surrogate slope contrasts at $2T$; only parameter grids passing every covariance check are eligible for mechanistic interpretation.')
    failed_notice = ('Every reported binary-control cell passes the retained validation tolerance.'
                     if not failed else 'Parameter grids containing a failed covariance check are diagnostic only and are excluded from mechanistic conclusions.')
    quarter = next(c for c in binary_rows if c['theta']==.25)
    qci = np.array(quarter['ci95'])
    moment = cells[1]
    mci = moment['four_sign_moment_difference_ci95']
    if not eligible[.25] or not qci[0,1] <= 0 <= qci[1,1] or mci[0] <= 0 <= mci[1]:
        raise ValueError('New results require revising the abstract/conclusion before replacing the PDF')
    refinement = json.loads((ROOT/'docs/innovation_20261031/refinement_report.json').read_text(encoding='utf-8'))
    refinement_text = '; '.join('$L='+str(r['spec']['L'])+'$: $'+
        f'{r["previous_gate"]["simultaneous_upper95"]:.5f}'+r'\to'+
        f'{r["refined_gate"]["simultaneous_upper95"]:.5f}'+'$' for r in refinement)
    extra = [c['binary_calibration']['additional_covariance_validation'] for c in result['cells']
             if 'additional_covariance_validation' in c.get('binary_calibration', {})]
    extra_notice = ''
    if extra:
        e = extra[0]
        extra_notice = (r'At $\theta=0.40$, $L=512$, the initial upper bound was '+
            f'{e["previous_gate"]["simultaneous_upper95"]:.5f}'+
            '. Keeping the fitted spectrum and all interface records fixed, we performed one '
            'separate covariance-only validation with 128 fresh runs and 64 sites per input '
            'and run. Its upper bound is '+f'{e["gate"]["simultaneous_upper95"]:.5f}'+
            '. These records were not used to fit the spectrum or interface slopes. '
            'No further validation retries or tolerance relaxation were used.')
    numbers = {
        'InnovationSizeRange': str(sizes[0])+r'\mbox{--}'+str(sizes[-1]),
        'InnovationSizeCount': str(len(sizes)),
        'InnovationCells': str(len(result['cells'])),
        'InnovationBinaryPass': str(len(result['cells'])-len(failed)),
        'InnovationStageNotice': stage,
        'InnovationFailedCells': failed_text or 'none',
        'InnovationRawResidualMin': f'{min(raw):.3f}',
        'InnovationRawResidualMax': f'{max(raw):.3f}',
        'InnovationSurrogateFinding': finding,
        'InnovationFailedInferenceNotice': failed_notice,
        'InnovationRefinementBounds': refinement_text,
        'InnovationLargeValidationNotice': extra_notice,
        'InnovationQuarterSlopeContrast': f'{quarter["delta_alpha"][1]:.4f}'+
            r'\;['+f'{qci[0,1]:.4f},{qci[1,1]:.4f}'+']',
        'InnovationQuarterMomentContrast': f'{moment["four_sign_moment_difference"]:.4f}'+
            r'\;['+f'{mci[0]:.4f},{mci[1]:.4f}'+']',
    }
    if failed:
        numbers['InnovationFailedInferenceNotice'] += ' Remaining failed cells: '+', '.join(
            r'$\theta='+str(c['theta'])+r'$, $L='+str(c['L'])+'$ (upper bound '+
            f'{c["binary_covariance_gate"]["simultaneous_upper95"]:.5f}'+')' for c in failed)+'.'
    lines = ['\\newcommand{\\'+key+'}{'+value+'}' for key, value in numbers.items()]
    rows = []
    for theta in [.1, .25, .4]:
        for key, label in [('gaussian_innovation', 'G independent'),
                           ('uniform_innovation', 'U independent'),
                           ('binary_match_vs_heldout_uniform', 'G binary target')]:
            c = next(c for c in comparisons if c['theta']==theta and c['comparison']==key)
            ci = np.asarray(c['ci95'])
            marker = r'^{\dagger}' if key=='binary_match_vs_heldout_uniform' and not eligible[theta] else ''
            rows.append('$'+f'{theta:.2f}'+marker+'$ & '+label+' & $'+f'{c["delta_alpha"][1]:.4f}'+
                        '$ & $['+f'{ci[0,1]:.4f},{ci[1,1]:.4f}'+r']$ \\')
    lines.append(r'\newcommand{\InnovationContrastRows}{'+'\n'+'\n'.join(rows)+'\n}')
    (ROOT/'manuscript/innovation_numbers.tex').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    (figdir/'sources.json').write_text(json.dumps(dict(summary=summary.relative_to(ROOT).as_posix(),
        summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest(), sizes=sizes, outputs=outputs), indent=2), encoding='utf-8')
    build_data_archive(summary, result, sizes)
    build = ROOT/'docs/innovation_20261031/pdf_build'
    build.mkdir(exist_ok=True)
    env = os.environ.copy()
    env['TECTONIC_CACHE_DIR'] = str(ROOT/'.dsh_tmp/tectonic/cache')
    command = [str(ROOT/'.dsh_tmp/tectonic/tectonic.exe'), '--only-cached',
               '--keep-logs', '--keep-intermediates', '--outdir', str(build), 'main.tex']
    run = subprocess.run(command, cwd=ROOT/'manuscript', env=env, capture_output=True,
                          encoding='utf-8', errors='replace')
    (build/'compiler_output.txt').write_text(run.stdout+run.stderr, encoding='utf-8')
    if run.returncode:
        raise RuntimeError('PDF compilation failed; see pdf_build/compiler_output.txt')
    import fitz
    with fitz.open(build/'main.pdf') as document:
        text = ' '.join(' '.join(page.get_text().split()) for page in document)
        if '??' in text or '[?]' in text or 'Stationary innovation and covariance-matching' not in text:
            raise RuntimeError('PDF content or cross-reference validation failed')
        if 'Overfull' in (build/'main.log').read_text(encoding='utf-8', errors='replace'):
            raise RuntimeError('PDF has overfull boxes; review the compile log')
        pdf_pages = len(document)
    target = ROOT/'manuscript/BD_noise_controls_updated.pdf'
    shutil.copyfile(build/'main.pdf', target)
    for name in ['main.pdf', 'BD_figures_revised.pdf']:
        try:
            shutil.copyfile(target, ROOT/'manuscript'/name)
        except PermissionError:
            print(f'Reader locked {name}; updated named PDF is available.', flush=True)
    subprocess.run([sys.executable, str(ROOT/'code/build_overleaf.py')], cwd=ROOT, check=True)
    build_submission_bundle()
    verification = dict(pdf=target.relative_to(ROOT).as_posix(), pages=pdf_pages,
                         sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                         control_sizes=sizes, unresolved_references=0, overfull_boxes=0,
                         summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest())
    (build/'verification.json').write_text(json.dumps(verification, indent=2), encoding='utf-8')
    print(f'UPDATED PDF: {target}; control sizes={sizes}', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summary', type=Path, required=True)
    ap.add_argument('--wait-for-full', action='store_true')
    ap.add_argument('--figure-only', action='store_true')
    args = ap.parse_args()
    if args.wait_for_full:
        ready = ROOT/'docs/innovation_20261031/results/completion.json'
        while not ready.exists():
            time.sleep(15)
    export(args.summary.resolve(), args.figure_only)


if __name__ == '__main__':
    main()
