"""Run bootstrap for exact binary, spectral residual and joint limit controls."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from controlled_analysis import adjacent_contrast
from strict_binary import OUT as RAW, ROOT

OUT = ROOT/'data/analysis/strict_binary_20261117'
SIZES = [64,128,256,512]
B = 4000


def report(estimate, draws, simultaneous=False):
    result = dict(estimate=np.asarray(estimate).tolist(),
                  ci95=np.quantile(draws,[.025,.975],axis=0).tolist())
    if simultaneous:
        sd = draws.std(axis=0)
        z = np.divide(draws-estimate,sd,out=np.zeros_like(draws),where=sd>0)
        critical = np.quantile(np.max(abs(z.reshape(B,-1)),axis=1),.95)
        result['simultaneous95'] = np.stack([estimate-critical*sd, estimate+critical*sd]).tolist()
    return result


def means(values, rng):
    draws = rng.multinomial(len(values),np.full(len(values),1/len(values)),size=B)/len(values)
    return values.mean(axis=0),np.tensordot(draws,values,axes=1)


def comparison(a,b,rng,square_root=True):
    x,xx = means(a,rng)
    y,yy = means(b,rng)
    if square_root:
        return np.sqrt(y/x),np.sqrt(yy/xx)
    return y-x,yy-xx


def main():
    rng = np.random.default_rng(2026111720)
    sources = {}
    def load(path):
        sources[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path,allow_pickle=False) as z:
            return dict(z)

    early, late, spectral, diagnostic, records = [],[],[],[],[]
    for sites in SIZES:
        z = load(RAW/f'L{sites}.npz')
        assert json.loads(str(z['spec']))['runs']==128
        records.append(z['late_w2'])
        early.append(comparison(*z['early_w2'].transpose(1,0,2),rng))
        late.append(comparison(*z['late_w2'].transpose(1,0,2),rng))
        # Three fixed Fourier bands summarise the low, intermediate and high modes.
        s=z['structure_factor'][:,:,1,:sites//2+1]
        bands=np.stack([s[:,:,lo:hi].sum(axis=-1) for lo,hi in [(1,5),(5,17),(17,sites//2+1)]],axis=-1)
        spectral.append(comparison(bands[:,0],bands[:,1],rng))
        d={}
        for key in ['deposition','corr','fourth','spatial']:
            d[key]=[report(*means(z[key][:,j],rng)) for j in [0,1]]
        diagnostic.append(d)
    def family(pairs, squared=False):
        y=np.stack([p[0] for p in pairs])
        b=np.stack([p[1] for p in pairs],axis=1)
        return report(y*y,b*b,True) if squared else report(y,b,True)

    alpha = adjacent_contrast(SIZES,[z[:,0] for z in records],[z[:,1] for z in records],draws=B)
    exact_path=RAW/'exact.json'
    sources[exact_path.relative_to(ROOT).as_posix()]=hashlib.sha256(exact_path.read_bytes()).hexdigest()
    result=dict(sizes=SIZES, exact=json.loads(exact_path.read_text()),
                early=family(early),late=family(late),alpha=alpha,spectral_bands=family(spectral,True),
                diagnostics=diagnostic)
    if '--parity-only' in sys.argv:
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/'parity.json').write_text(json.dumps(result,indent=2)+'\n')
        return

    residual, residual_sizes = [],[64,128,256]
    residual_w, residual_draws = [],[]
    for L in residual_sizes:
        z=load(RAW/'residual'/f'L{L}.npz')
        w,wb=means(z['late_w2'],rng)
        log,logb=.5*np.log(w),.5*np.log(wb)
        residual_w.append(log);residual_draws.append(logb)
        cov,covb=means(z['covariance'][:,1:]-z['covariance'][:,:1],rng)
        pred=z['theoretical_sign'][1:]-z['theoretical_sign'][:1]
        lags=z['lags']
        residual.append(dict(L=L,width=report(np.exp(log[1:]-log[:1]),np.exp(logb[:,1:]-logb[:,:1])),
                             all_lag_displacement=np.max(abs(pred),axis=1).tolist(),
                             sampled_lags=lags.tolist(),predicted=pred[:,lags].tolist(),
                             measured=report(cov,covb),
                             max_measured_prediction_gap=np.max(abs(cov-pred[:,lags]),axis=1).tolist(),
                             paired_K4=report(*means(z['four_cumulant'][:,1:]-z['four_cumulant'][:,:1],rng))))
    log=np.stack(residual_w)
    boot=np.stack(residual_draws,axis=1)
    da=np.diff(log[:,1:]-log[:,:1],axis=0)/np.log(2)
    dab=np.diff(boot[:,:,1:]-boot[:,:,:1],axis=1)/np.log(2)
    result['residual']=residual
    result['residual_alpha']=report(da,dab,True)

    joint_groups, joint = [],[]
    for j in range(3):
        cells=[]
        for L in [64,128]:
            z=load(RAW/f'joint{j}'/f'theta0.25_L{L}.npz')
            spec=json.loads(str(z['spec']))
            assert spec['runs']==32
            w,wb=means(z['late_w2'],rng)
            cells.append((w,wb))
            joint.append(dict(level=j,L=L,spec=spec,
                              width=report(np.sqrt(w),np.sqrt(wb)),
                              latest_duration_ratio=report(np.sqrt(w[:,-1]/w[:,-2]),np.sqrt(wb[:,:,-1]/wb[:,:,-2]))))
        (a,ab),(b,bb)=cells
        joint_groups.append(cells)
    for j in range(1,3):
        for i,L in enumerate([64,128]):
            a,ab=joint_groups[j-1][i]
            b,bb=joint_groups[j][i]
            row=next(x for x in joint if x['L']==L and x['level']==j)
            row['latest_vs_previous']=report(np.sqrt(b[:,-1]/a[:,-1]),np.sqrt(bb[:,:,-1]/ab[:,:,-1]))
            row['matched_vs_previous']=report(np.sqrt(b[:,-2]/a[:,-1]),np.sqrt(bb[:,:,-2]/ab[:,:,-1]))
    local=[(.5*np.log(b[0][:,-1]/a[0][:,-1])/np.log(2), .5*np.log(b[1][:,:,-1]/a[1][:,:,-1])/np.log(2)) for a,b in joint_groups]
    result['joint']=joint
    result['joint_alpha']=[report(*p) for p in local]
    result['joint_alpha_change']=report(local[2][0]-local[1][0],local[2][1]-local[1][1])
    result['bootstrap_draws']=B
    result['bootstrap_seed']=2026111720
    result['sources']=sources
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Analysed strict binary, residual and joint cohorts.',flush=True)


if __name__=='__main__':
    main()
