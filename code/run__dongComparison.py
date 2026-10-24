"""Matched phase, section vorticity, circulation and pressure-jump comparisons."""
from pathlib import Path
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import save_hybrid


def pressure_jump(data,j):
    xy=data[f'{j}_actual_outline'];cp=data[f'{j}_actual_pstar']
    edge=np.diff(xy,axis=0);mid=(xy[1:]+xy[:-1])/2
    orientation=np.sign(np.sum(xy[:-1,0]*xy[1:,1]-xy[:-1,1]*xy[1:,0]))
    upper=-orientation*edge[:,0]>0
    assert np.mean(mid[upper,1]) > np.mean(mid[~upper,1])
    values=cp[:-1]+cp[1:]
    query=np.linspace(xy[:,0].min()+.05*np.ptp(xy[:,0]),xy[:,0].max()-.05*np.ptp(xy[:,0]),91)
    profiles=[]
    for mask in [upper,~upper]:
        order=np.argsort(mid[mask,0]);profiles.append(np.interp(query,mid[mask,0][order],values[mask][order]))
    return query,profiles[1]-profiles[0]


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    reports=Path(__file__).resolve().parents[1]/'reports'
    out=reports/'2026-10-23-dong-comparison';out.mkdir(parents=True,exist_ok=True)
    u=.29;span=np.array([.2,.3,.4,.5,.7,.9]);layers=1
    signals=np.load(reports/'2026-10-23-lift-cause-controls/cause_signals.npz')
    force=signals['dynamic_contrast'];peaks=[int(np.argmax(force[:8])),8+int(np.argmax(force[8:16]))]
    phases=sorted(set((p+d)%16 for p in peaks for d in [-1,0,1]))
    records=[];manifest=[]
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    for phase in phases:
        caches={};jsons={}
        for case in ['baseline','isolated']:
            path=study/'vortex_structures_v4'/f'{case}_{phase:02d}_sections.npz'
            caches[case]=dict(np.load(path));manifest.append(dict(name=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            jsons[case]=json.loads((path.with_name(f'{case}_{phase:02d}.json')).read_text(encoding='utf-8'))
        row=dict(phase=phase/16,is_force_peak=phase in peaks,paired_two_fin_increment=float(force[phase]),sections=[])
        for j,s in enumerate(span):
            entry=dict(span=float(s),cases={})
            key=f'{j}_{layers}_';b=caches['baseline'];i=caches['isolated']
            for name in ['x','y','upper','lower','basis']:
                np.testing.assert_allclose(b[key+name],i[key+name],atol=1e-10)
            xx,yy=np.meshgrid(b[key+'x'],b[key+'y'])
            xc=(xx-b[key+'leading_x'])/b[key+'chord']
            gap=yy-b[key+'upper'][None,:]
            window=(xc>=0)&(xc<=.5)&(gap>=.012)&(gap<=.04)
            common=window & b[key+'finite'] & i[key+'finite']
            entry['fixed_window_coverage']=float(common.sum()/max(window.sum(),1))
            for case,data in caches.items():
                query,jump=pressure_jump(data,j)
                shifted=dict(data);shifted[f'{j}_actual_pstar']=data[f'{j}_actual_pstar']+.73
                np.testing.assert_allclose(pressure_jump(shifted,j)[1],jump,atol=1e-12)
                candidates=[r for r in jsons[case]['records'] if np.isclose(r['span'],s) and r['side']=='upper' and r['zone']=='front' and r['excluded_layers']==layers]
                entries=[dict(threshold=r['threshold'],gamma_over_U=r.get('circulation_over_UL'),closed=r.get('closed_resolved_contour',False),found=r['found']) for r in candidates]
                area=float(data[key+'spacing'])**2
                entry['cases'][case]=dict(mean_jump=float(np.trapezoid(jump,query)/np.ptp(query)),circulation_candidates=entries,
                    fixed_window_signed_gamma_over_U=float(np.sum(data[key+'omega_span'][common])*area/u),
                    fixed_window_negative_gamma_over_U=float(np.sum(np.minimum(data[key+'omega_span'][common],0))*area/u),
                    fixed_window_mean_span_velocity_over_U=float(np.mean(data[key+'velocity_local'][...,2][common])/u))
            entry['delta_jump']=entry['cases']['baseline']['mean_jump']-entry['cases']['isolated']['mean_jump']
            records.append(dict(phase=phase/16,**entry));row['sections'].append(entry)
        if phase not in peaks:continue
        fig,axs=plt.subplots(2,3,figsize=(10,5.6),layout='constrained')
        vmax=max(float(np.nanquantile(abs(caches[c][f'{j}_{layers}_omega_span'])/u,.99)) for c in caches for j in [0,2,4])
        for r,(case,data) in enumerate(caches.items()):
            for col,j in enumerate([0,2,4]):
                key=f'{j}_{layers}_';ax=axs[r,col];chord=float(data[key+'chord']);le=float(data[key+'leading_x'])
                x=(data[key+'x']-le)/chord;y=data[key+'y']/chord
                im=ax.pcolormesh(x,y,data[key+'omega_span']/u,cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='auto',zorder=0)
                vel=data[key+'velocity_local'];xx,yy=np.meshgrid(x,y)
                ax.quiver(xx[::7,::7],yy[::7,::7],vel[::7,::7,0]/u,vel[::7,::7,1]/u,scale=20,width=.004,zorder=0)
                poly=data[key+'outline'];ax.fill((poly[:,0]-le)/chord,poly[:,1]/chord,color='.65',zorder=0)
                ax.set(title=f'{case}, span={span[j]:.1f}',xlabel='Chord coordinate x/c',ylabel='Normal coordinate y/c',xlim=(-.2,1.3))
                ax.set_rasterization_zorder(1)
        fig.suptitle(f'Matched sections at t/T={phase/16:.4f}; local velocity and spanwise vorticity')
        fig.colorbar(im,ax=axs,label='Local spanwise vorticity / U',shrink=.8)
        save_hybrid(fig,out/f'18_sections_phase_{phase:02d}',600)
        fig,axs=plt.subplots(1,3,figsize=(10,3.2),layout='constrained')
        for case,color,label in [('baseline','#D55E00','Body + fin'),('isolated','#0072B2','Only fin')]:
            entries=[s['cases'][case] for s in row['sections']]
            axs[0].plot(span,[s['mean_jump'] for s in entries],'o-',color=color,label=label,zorder=0)
            axs[2].plot(span,[-v['fixed_window_negative_gamma_over_U'] for v in entries],'o-',color=color,label=label,zorder=0)
        axs[1].plot(span,[s['delta_jump'] for s in row['sections']],'o-',color='black',zorder=0)
        for ax in axs:
            ax.set(xlabel='Span fraction');ax.set_rasterization_zorder(1);ax.spines[['top','right']].set_visible(False)
        axs[0].set(ylabel='Chord-mean lower-minus-upper Cp');axs[0].legend(frameon=False)
        axs[1].set(ylabel='Difference in pressure jump, B - I')
        axs[2].set(ylabel='Negative vorticity integral magnitude / U',title='Matched upper-front fluid window')
        fig.suptitle(f'Pressure and circulation compared at t/T={phase/16:.4f}')
        save_hybrid(fig,out/f'19_pressure_circulation_phase_{phase:02d}',600)
    for svg in out.glob('*.svg'):
        svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    result=dict(peaks=[p/16 for p in peaks],selection='Maximum material-paired two-fin pressure increment in each half cycle; neighboring saved phases included.',records=records,provenance=manifest,
        caution='Candidate front-upper signed-vorticity connected components, not automatically isolated LEVs. Truncated contours are retained and flagged; do not equate their circulation with true total vortex circulation. Geometry-aligned slices, chord mean pressure over 5–95%.')
    (out/'dong_comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(peaks=result['peaks'],peak_records=[r for r in records if r['phase'] in result['peaks']]),indent=2))


if __name__=='__main__':main()
