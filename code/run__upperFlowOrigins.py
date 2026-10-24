"""Trace fluid in the newly measured upper-front enhancement region backwards."""
from pathlib import Path
from functools import lru_cache
import sys
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    archive=study/'code_archive/20261023_before_cleanup'
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-dong-evolution'
    sys.path.insert(0,str(archive))
    from trace__temporal_source_paths import trajectories
    from read_volume import nearest_surface
    cases=['baseline','isolated'];phase=12;rows=[]
    reference=dict(np.load(study/'vortex_structures_v4'/f'baseline_{phase:02d}_sections.npz'))
    seeds=[];labels=[]
    for j in [0,1,2]:
        key=f'{j}_1_';xx,yy=np.meshgrid(reference[key+'x'],reference[key+'y'])
        xc=(xx-reference[key+'leading_x'])/reference[key+'chord'];gap=yy-reference[key+'upper'][None,:]
        eligible=(xc>=0)&(xc<=.5)&(gap>=.012)&(gap<=.04)&reference[key+'finite']&(reference[key+'Q']>0)
        strength=np.where(eligible,-reference[key+'omega_span'],-np.inf)
        iy,ix=np.unravel_index(np.argmax(strength),strength.shape)
        seed=reference[key+'origin']+reference[key+'basis']@np.array([xx[iy,ix],yy[iy,ix],0.])
        for offset in [[0,0,0],[.0015,0,0],[-.0015,0,0],[0,.0015,0],[0,-.0015,0],[0,0,.0015],[0,0,-.0015]]:
            seeds.append(seed+offset);labels.append(dict(span=[.2,.3,.4][j],offset=offset))
    for case in cases:
        @lru_cache(maxsize=3)
        def frame(index):
            with np.load(study/'vortex_source_v8/volumes'/f'{case}_{index:02d}.npz') as d:
                axes=(d['z'],d['y'],d['x'])
                return {k:RegularGridInterpolator(axes,d[field].astype('f4'),bounds_error=False,fill_value=(0 if k=='mask' else np.nan)) for k,field in [('velocity','velocity'),('mask','fluid'),('Q','Q'),('omega','omega')]}
        def series(key):return [lambda points,index=j:frame(index)[key](points) for j in range(17)]
        geometry=dict(np.load(study/'loads'/f'{case}.npz'))
        for steps in [512,1024]:
            histories,reasons=trajectories(np.array(seeds),series('velocity'),series('mask'),series('Q'),series('omega'),steps,start=phase/16)
            for label,history,reason in zip(labels,histories,reasons):
                endpoint=history[-1];t=endpoint['t'];j=min(int(t*16),15);a=t*16-j
                points=(1-a)*geometry['points_history'][j]+a*geometry['points_history'][j+1]
                ids,bary,distance=nearest_surface(points,geometry['faces'],np.array(endpoint['xyz'])[None])
                rows.append(dict(case=case,steps=steps,**label,reason=reason,history=history,distance=float(distance[0]),nearest_fin=bool(geometry['pectoral_mask'][ids[0]])))
            print(case,steps,'completed',flush=True)
        frame.cache_clear()
    out.mkdir(parents=True,exist_ok=True)
    summary=[]
    for case in cases:
        selected=[r for r in rows if r['case']==case and r['steps']==1024]
        coarse=[next(q for q in rows if q['case']==case and q['steps']==512 and q['span']==r['span'] and q['offset']==r['offset']) for r in selected]
        summary.append(dict(case=case,total=len(selected),within_006_of_fin=sum(r['distance']<.006 and r['nearest_fin'] for r in selected),
            within_006_of_body=sum(r['distance']<.006 and not r['nearest_fin'] for r in selected),
            max_endpoint_difference_between_step_sizes=max(float(np.linalg.norm(np.array(r['history'][-1]['xyz'])-q['history'][-1]['xyz'])) for r,q in zip(selected,coarse)),
            max_stop_time_difference=max(abs(r['history'][-1]['t']-q['history'][-1]['t']) for r,q in zip(selected,coarse))))
    (out/'upper_flow_origins.json').write_text(json.dumps(dict(records=rows,summary=summary,definition='Identical paired release positions at t/T=.75; baseline Q-positive maximum within each common upper-front window plus six small offsets. Backward fluid paths, not vortex identities; near-wall proximity is not a demonstrated vorticity source.'),indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from fn___plotEvidence import save_hybrid
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(8,3.4),layout='constrained')
    for ax,vertical in zip(axes,[1,2]):
        for case,color in [('baseline','#D55E00'),('isolated','#0072B2')]:
            for r in rows:
                if r['case']!=case or r['steps']!=1024 or r['offset']!=[0,0,0]:continue
                points=np.array([p['xyz'] for p in r['history']])
                ax.plot(points[:,0],points[:,vertical],color=color,label=case if r['span']==.2 else '_nolegend_',zorder=0)
                ax.scatter(points[0,0],points[0,vertical],marker='s',color=color,s=20,zorder=0)
                ax.scatter(points[-1,0],points[-1,vertical],marker='x',color=color,s=25,zorder=0)
        ax.set(xlabel='x/L',ylabel=['','y/L','z/L'][vertical]);ax.set_rasterization_zorder(1)
        ax.spines[['top','right']].set_visible(False)
    axes[0].legend(frameon=False)
    fig.suptitle('Backward fluid paths: squares = release at 0.75T; crosses = stopping points')
    save_hybrid(fig,out/'22_upper_flow_paths',600)
    svg=out/'22_upper_flow_paths.svg';svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(json.dumps([dict(case=r['case'],span=r['span'],stop_phase=r['history'][-1]['t'],distance=r['distance'],nearest_fin=r['nearest_fin'],reason=r['reason']) for r in rows if r['steps']==1024 and r['offset']==[0,0,0]],indent=2))


if __name__=='__main__':main()
