"""Plot the strictly pairwise-matched growth experiment from saved runs."""
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from refocus_figures import COLORS, WIDTH, EvidenceFigures, error, panel
from strict_binary_analysis import OUT, RAW, SIZES, means


def main():
    root = RAW.parents[2]
    figs = EvidenceFigures(root/'figs/strict_binary_20261117',root/'docs/strict_binary_20261117')
    summary = figs.read_json((OUT/'summary.json').relative_to(root).as_posix())
    for rel, sha in summary['sources'].items():
        figs.register(rel, sha)
    z = figs.read_npz((RAW/'L512.npz').relative_to(root).as_posix())
    fig, axes = plt.subplots(2,2,figsize=(WIDTH,5.5))
    fig.subplots_adjust(left=.105,right=.98,top=.91,bottom=.16,hspace=.64,wspace=.31)
    a,b,c,d = axes.flat
    name='fig4_strict_binary'
    rng=np.random.default_rng(2026111721)
    for i,key,label in [(0,'G','Independent'),(1,'U','Parity block')]:
        time=z['w_times'][0,i,::4]
        value,boot=means(z['w_series'][:,i,::4],rng)
        band=np.quantile(np.sqrt(boot),[.025,.975],axis=0)
        a.plot(time,np.sqrt(value),color=COLORS[key],label=label)
        a.fill_between(time,*band,color=COLORS[key],alpha=.15,lw=0)
        figs.add(name,'a',label,time,np.sqrt(value),band,theta=0,size=512,n=128)
    a.set(xscale='log',yscale='log',xlabel='Observation time, t',ylabel=r'RMS width, $W(t)$')
    a.legend(loc='lower right',fontsize=8)
    panel(a,'a','Width evolution, L = 512')
    late=summary['late']
    for i,key,label in [(0,'G','T'),(1,'U','2T')]:
        x=np.arange(4)+(i-.5)*.17
        y=np.array(late['estimate'])[:,i]
        bounds=np.array(late['simultaneous95'])[:,:,i]
        error(b,x,y,bounds,COLORS[key],label,['o','s'][i],ls='none')
        figs.add(name,'b',label,SIZES,y,bounds,theta=0,n=128)
    b.set(xticks=range(4),xticklabels=SIZES,xlabel='System size, L',ylabel=r'$W_B/W_A$')
    b.axhline(1,color='.5',lw=.7)
    b.legend(loc='lower right',ncol=2)
    panel(b,'b','Late-window width ratio')
    alpha=summary['alpha']
    for i,key,label in [(0,'G','T'),(1,'U','2T')]:
        x=np.arange(3)+(i-.5)*.16
        y=np.array(alpha['delta'])[:,i]
        bounds=np.array(alpha['ci95'])[:,:,i]
        error(c,x,y,bounds,COLORS[key],label,['o','s'][i],ls='none')
        figs.add(name,'c',label,SIZES[:-1],y,bounds,theta=0,n=128)
    c.set(xticks=range(3),xticklabels=['64→128','128→256','256→512'],
          xlabel='Adjacent sizes',ylabel=r'$\alpha_B-\alpha_A$')
    c.axhline(0,color='.5',lw=.7)
    panel(c,'c','Effective-exponent contrast')
    for i,key,label in [(0,'G','Independent'),(1,'U','Parity block')]:
        q=2*np.pi*np.arange(1,257)/512
        value,boot=means(z['structure_factor'][:,i,1,1:257],rng)
        bounds=np.quantile(boot,[.025,.975],axis=0)
        assert np.min(bounds)>0
        d.plot(q,value,color=COLORS[key],label=label)
        d.fill_between(q,*bounds,color=COLORS[key],alpha=.18,lw=0)
        figs.add(name,'d',label,q,value,bounds,theta=0,size=512,n=128)
    d.set(xscale='log',yscale='log',xlabel='Wave number, q',ylabel=r'$S(q)$')
    panel(d,'d','Spatial spectrum, L = 512')
    fig.text(.105,.072,'128 independent runs per input and size. (b) Simultaneous 95% intervals over eight width ratios.',fontsize=8)
    fig.text(.105,.034,'Other bars and bands: pointwise 95% whole-run bootstrap intervals. (d) Snapshots at 2T.',fontsize=8)
    figs.save(fig,name)
    print('Created Fig. 4 and panel source data.')


if __name__=='__main__':
    main()
