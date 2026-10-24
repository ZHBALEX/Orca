"""Measure rootward fluid flux below the isolated-fin reference geometry."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import save_hybrid


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    gates=[5.895,5.905,5.915];u=.29;rows=[]
    for phase in range(17):
        data={case:dict(np.load(study/'vortex_source_v8/volumes'/f'{case}_{phase:02d}.npz')) for case in ['baseline','isolated']}
        b,i=data.values()
        x,y=i['x'],i['y'];xx,yy=np.meshgrid(x,y);area=np.outer(np.gradient(y),np.gradient(x))
        bmask=RegularGridInterpolator((b['z'],b['y'],b['x']),b['fluid'].astype(float),method='nearest',bounds_error=True)
        bw=RegularGridInterpolator((b['z'],b['y'],b['x']),np.where(b['fluid'],b['velocity'][...,2],0)/u,bounds_error=True)
        for target in gates:
            iz=int(np.argmin(abs(i['z']-target)))
            solid=~i['fluid'][iz].astype(bool)
            selected=solid & (xx>10.15)&(xx<10.32)&(yy>5.88)&(yy<6.01)
            has=selected.any(0);lower=np.min(np.where(selected,yy,np.inf),axis=0)
            for h0,h1 in [(.005,.02),(.008,.025)]:
                gate=has[None,:]&(yy<lower[None,:]-h0)&(yy>lower[None,:]-h1)
                assert gate.any(), (phase,target,h0)
                record=dict(phase=phase/16,target_z=target,z=float(i['z'][iz]),h0=h0,h1=h1,area=float(area[gate].sum()),cases={})
                for case,d in data.items():
                    if case=='baseline':
                        query=np.stack([np.full(gate.sum(),i['z'][iz]),yy[gate],xx[gate]],axis=-1)
                        fluid_here=bmask(query)>.5;w_here=bw(query);weights=area[gate]
                    else:
                        fluid_here=d['fluid'][iz][gate].astype(bool);w_here=d['velocity'][iz,...,2][gate]/u;weights=area[gate]
                    w=w_here[fluid_here];weights=weights[fluid_here]
                    assert np.isfinite(w).all()
                    record['cases'][case]=dict(rootward_flux=float(np.sum(np.maximum(w,0)*weights)),outward_flux=float(np.sum(np.minimum(w,0)*weights)),
                        net_flux=float(np.sum(w*weights)),fluid_gate_fraction=float(weights.sum()/area[gate].sum()))
                rows.append(record)
    summaries=[]
    for target in gates:
        for h0,h1 in [(.005,.02),(.008,.025)]:
            selected=[r for r in rows if r['target_z']==target and r['h0']==h0];t=np.array([r['phase'] for r in selected]);s=dict(target_z=target,h0=h0,h1=h1,cases={})
            for case in ['baseline','isolated']:
                s['cases'][case]={key:float(np.trapezoid([r['cases'][case][key] for r in selected],t)) for key in ['rootward_flux','outward_flux','net_flux','fluid_gate_fraction']}
            s['rootward_flux_ratio']=s['cases']['baseline']['rootward_flux']/s['cases']['isolated']['rootward_flux'];summaries.append(s)
    out.mkdir(parents=True,exist_ok=True)
    (out/'root_gate_flux.json').write_text(json.dumps(dict(rows=rows,summary=summaries,definition='Integral of W/U over moving-y lower-side windows at fixed Cartesian z. Positive W is toward the root/body. Gate bounds follow the isolated-fin voxel surface at each phase. Solid portions carry no fluid flux; this is an open gate, not a closed control-volume budget or a count of completed lower-to-upper paths.'),indent=2),encoding='utf-8')
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(8,3.2),layout='constrained')
    selected=[r for r in rows if r['target_z']==5.915 and r['h0']==.005];t=[r['phase'] for r in selected]
    for case,color in [('baseline','#D55E00'),('isolated','#0072B2')]:
        for ax,key in zip(axes,['rootward_flux','net_flux']):
            ax.plot(t,[r['cases'][case][key]*1e4 for r in selected],'o-',color=color,label=case,zorder=0,ms=3)
    for ax,title in zip(axes,['Rootward positive flux','Signed lateral flux']):
        ax.set(xlabel='t/T',ylabel='Flux / (U L²) × 10000',title=title);ax.set_rasterization_zorder(1);ax.axhline(0,color='.65',lw=.6,zorder=0)
    axes[0].legend(frameon=False)
    save_hybrid(fig,out/'25_root_gate_flux',600)
    svg=out/'25_root_gate_flux.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
