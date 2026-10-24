"""Body/fin crossflow sequence, retaining Cartesian directions and solid masks."""
from pathlib import Path
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import save_hybrid


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow';out.mkdir(parents=True,exist_ok=True)
    phases=[6,8,10,12];u=.29;position=10.24
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(2,4,figsize=(12,5),layout='constrained')
    for row,case in enumerate(['baseline','isolated']):
        for col,phase in enumerate(phases):
            with np.load(study/'vortex_source_v8/volumes'/f'{case}_{phase:02d}.npz') as d:
                ix=int(np.argmin(abs(d['x']-position)));z,y=d['z'],d['y'];v=d['velocity'][:,:,ix];mask=d['fluid'][:,:,ix].astype(bool)
                zz,yy=np.meshgrid(z,y,indexing='ij');ax=axes[row,col]
                im=ax.pcolormesh(z,y,np.where(mask,v[...,2]/u,np.nan).T,cmap='RdBu_r',vmin=-.35,vmax=.35,shading='auto',zorder=0)
                ax.contourf(z,y,mask.T,levels=[-.5,.5],colors=['.5'],zorder=0)
                ax.quiver(zz[::3,::3],yy[::3,::3],np.where(mask,v[...,2]/u,np.nan)[::3,::3],np.where(mask,v[...,1]/u,np.nan)[::3,::3],scale=6,width=.005,zorder=0)
                ax.set(xlim=(5.85,5.96),ylim=(5.88,6.015),xlabel='z/L (toward body →)',ylabel='y/L',title=f'{case}; t/T={phase/16:g}')
                ax.set_rasterization_zorder(1)
    fig.colorbar(im,ax=axes,label='Cartesian lateral velocity W/U',shrink=.8)
    fig.suptitle(f'Body–fin crossflow near x/L={position}; grey = solid; arrows = (W,V)/U')
    save_hybrid(fig,out/'24_root_crossflow_sequence',600)
    rows=list(csv.DictReader((study/'section_metrics.csv').open(encoding='utf-8')))
    records=[]
    for span in [.2,.3,.4,.5,.7,.9]:
        result=dict(span=span,cases={})
        for case in ['baseline','isolated']:
            rr=sorted([r for r in rows if r['case']==case and np.isclose(float(r['span']),span)],key=lambda r:float(r['phase']))
            t=np.array([float(r['phase']) for r in rr])
            result['cases'][case]={key:dict(mean=float(np.trapezoid(np.array([float(r[key]) for r in rr]),t)),values=[float(r[key]) for r in rr]) for key in ['W_upper_20','W_lower_20','U_upper_20','U_lower_20','V_upper_20','V_lower_20']}
        records.append(result)
    decomposition=[]
    for phase in range(17):
        data={case:dict(np.load(study/'vortex_structures_v4'/f'{case}_{phase:02d}_sections.npz')) for case in ['baseline','isolated']}
        for j,span in enumerate([.2,.3,.4]):
            key=f'{j}_1_';b,i=data.values();xx,yy=np.meshgrid(b[key+'x'],b[key+'y']);xc=(xx-b[key+'leading_x'])/b[key+'chord'];gap=yy-b[key+'upper'][None,:]
            mask=(xc>=0)&(xc<=.5)&(gap>=.012)&(gap<=.04)&b[key+'finite']&i[key+'finite']
            delta=(b[key+'velocity_local']-i[key+'velocity_local'])@b[key+'basis'].T
            components=np.mean(delta[mask],axis=0)*b[key+'basis'][:,2]/u
            direct=np.mean((b[key+'velocity_local']-i[key+'velocity_local'])[...,2][mask])/u
            np.testing.assert_allclose(components.sum(),direct,atol=1e-12)
            decomposition.append(dict(phase=phase/16,span=span,local_span_direction=b[key+'basis'][:,2].tolist(),cartesian_velocity_difference=(delta[mask].mean(0)/u).tolist(),projection_components=components.tolist()))
    (out/'root_crossflow.json').write_text(json.dumps(dict(records=records,local_span_projection_audit=decomposition,
        scope='Cartesian velocity near fixed-Z sections and local-span projection at geometry-aligned sections are different samplings. Positive W on the negative-z fin points toward the body. These are velocity diagnostics, not integrated cross-root flux.'),indent=2),encoding='utf-8')
    for svg in out.glob('*.svg'):svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(json.dumps([dict(span=r['span'],**{case:{k:v['mean'] for k,v in r['cases'][case].items()} for case in r['cases']}) for r in records],indent=2))


if __name__=='__main__':main()
