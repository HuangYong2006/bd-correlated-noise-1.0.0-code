"""Paired filter inference and fixed-record duration controls for the PRE revision."""
from pathlib import Path
import hashlib
import json
import numpy as np
from bd_core import generate_spectral_eta, _simulate_from_field
from revision_figures import configure_style, canvas, save

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/analysis/submission_20261030'


def paired_slopes():
    rows, sources = [], []
    rng = np.random.default_rng(20261030)
    for theta in [0., .1, 1/6, .2, .25, .3, .4]:
        logs, draws = [], []
        for L in [32, 64, 128, 256, 512]:
            values, specs = [], []
            for tag in ['powerlaw', 'fgn']:
                # Archive filenames encode only the parameter's decimal point.
                p = ROOT / 'data/raw' / f'rev20261026_{tag}' / f'bd_p{theta:.4f}_L{L}.npz'.replace(f'{theta:.4f}', f'{theta:.4f}'.replace('.', 'd'))
                with np.load(p, allow_pickle=False) as z:
                    specs.append(json.loads(str(z['spec'])))
                    assert str(z['noise_layout']) == 'spectral-run-time-site-v2'
                    w = z['w_series'].astype(float)
                    values.append(w[:, -max(1, int(w.shape[1]*.2)):].mean(axis=1))
                sources.append({'file': p.relative_to(ROOT).as_posix(), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
            assert all(specs[0][k] == specs[1][k] for k in ['seed', 'L', 'T', 'n_runs'])
            idx = rng.integers(len(values[0]), size=(2000, len(values[0])))
            logs.append(.5*np.log(values[1].mean()/values[0].mean()))
            draws.append(.5*np.log(values[1][idx].mean(axis=1)/values[0][idx].mean(axis=1)))
        x = np.log([32,64,128,256,512]); x -= x.mean()
        delta = np.asarray(logs) @ x / (x @ x)
        boot = np.asarray(draws).T @ x / (x @ x)
        rows.append({'theta':theta, 'delta_alpha':float(delta), 'se':float(boot.std(ddof=1)),
                     'ci95':np.quantile(boot,[.025,.975]).tolist()})
    return rows, sources


def duration_controls():
    rows = []
    for theta in [.1, .4]:
        for L in [128,256]:
            path = OUT / f'fixed_record_theta{theta}_L{L}.npz'
            if not path.exists():
                old = ROOT / 'data/raw/rev20261026_fgn' / f'bd_p{theta:.4f}_L{L}.npz'.replace(f'{theta:.4f}', f'{theta:.4f}'.replace('.', 'd'))
                with np.load(old) as z: T = json.loads(str(z['spec']))['T']
                seed = 20261030 + L + int(theta*1000)
                all_widths = np.empty((3,32))
                for r0, eta in generate_spectral_eta(theta,L,4*T,32,seed,spectrum='fgn'):
                    for j, factor in enumerate([1,2,4]):
                        t = factor*T; stride = max(1,t//5000)
                        w = np.zeros((len(eta),t//stride),np.float32)
                        profiles = np.zeros((len(eta),1,L),np.int32)
                        _simulate_from_field(eta[:,:t,:],0,L,t,stride,np.array([t]),profiles,w)
                        assert np.max(np.abs(profiles-np.roll(profiles,1,axis=2))) <= 1
                        all_widths[j,r0:r0+len(eta)] = w[:,-max(1,int(w.shape[1]*.2)):].mean(axis=1)
                np.savez_compressed(path, squared_width=all_widths,T=T,L=L,theta=theta,
                                    seed=seed,n_runs=32,record_length=4*T)
            with np.load(path) as z:
                values=z['squared_width']; T=int(z['T'])
            rng=np.random.default_rng(20261030+L)
            idx=rng.integers(32,size=(2000,32))
            for j,factor in enumerate([2,4],1):
                b=np.sqrt(values[j][idx].mean(axis=1)/values[0][idx].mean(axis=1))
                rows.append({'theta':theta,'L':L,'T':T,'factor':factor,'n_runs':32,
                             'ratio':float(np.sqrt(values[j].mean()/values[0].mean())),
                             'se':float(b.std(ddof=1)),'ci95':np.quantile(b,[.025,.975]).tolist(),
                             'raw_file':path.relative_to(ROOT).as_posix(),
                             'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            print(f'Completed fixed-record control theta={theta} L={L}',flush=True)
    return rows


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    paired,sources=paired_slopes()
    (OUT/'paired.json').write_text(json.dumps({'rows':paired,'sources':sources},indent=2))
    duration=duration_controls()
    result={'paired':paired,'duration':duration,'bootstrap':2000,'seed':20261030,
            'duration_protocol':'One length-4T fGn-spectrum record; same run prefixes; flat initial interfaces; last 20% regularly sampled widths; fixed low-frequency cutoff.'}
    (OUT/'summary.json').write_text(json.dumps(result,indent=2))
    configure_style();fig,axes=canvas(height=80)
    axes[0].errorbar([r['theta'] for r in paired],[r['delta_alpha'] for r in paired],
                     yerr=[r['se'] for r in paired],fmt='o-',color='#1f6fb4',capsize=2,ms=3)
    axes[0].axhline(0,color='#4d4d4d',lw=.7)
    axes[0].set(xlabel='θ',ylabel=r'$\alpha_{\mathrm{fGn}}-\alpha_{\mathrm{power}}$',title='Paired filter comparison')
    for theta,color in [(.1,'#1f6fb4'),(.4,'#d1495b')]:
        for L,marker in [(128,'o'),(256,'s')]:
            data=[r for r in duration if r['theta']==theta and r['L']==L]
            axes[1].errorbar([r['factor'] for r in data],[r['ratio'] for r in data],
                             yerr=[r['se'] for r in data],color=color,marker=marker,ms=3,capsize=2,
                             label=f'θ={theta}, L={L}')
    axes[1].axhline(1,color='#4d4d4d',lw=.7)
    axes[1].set(xlabel='Observation duration / T',ylabel='Late-time width / width at T',title='Fixed length-4T input record',xticks=[2,4])
    axes[1].legend(fontsize=6,loc='best')
    save(fig,ROOT/'figs/revision_20261026','fig5_controls')
    lines=[]
    for r in paired:
        label='1/6' if abs(r['theta']-1/6)<1e-6 else f"{r['theta']:.2f}"
        lines.append(f"${label}$ & ${r['delta_alpha']:.5f}$ & $[{r['ci95'][0]:.5f},{r['ci95'][1]:.5f}]$ \\\\")
    (ROOT/'manuscript/submission_results.tex').write_text('\\newcommand{\\PairedRows}{\n'+'\n'.join(lines)+'\n}\n')


if __name__=='__main__':
    main()
