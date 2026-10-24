"""Matched moving/static contrasts for pressure and actual upper-fin windows."""
from pathlib import Path
import json
import csv
import numpy as np
from scipy.spatial import cKDTree
from run__dongComparison import pressure_jump


def top_at(x,outline):
    a,b=outline[:-1],outline[1:];dx=b[:,0]-a[:,0]
    frac=np.divide(x[...,None]-a[:,0],dx,out=np.full(x.shape+(len(a),),np.nan),where=abs(dx)>1e-12)
    yy=a[:,1]+frac*(b[:,1]-a[:,1]);valid=(frac>=0)&(frac<=1)
    return np.max(np.where(valid,yy,-np.inf),axis=-1)


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    spans=[.2,.3,.4,.5,.7,.9];u=.29;rows=[]
    for label,cases,frames in [('moving',['baseline','isolated'],range(17)),('static',['static_body','static_isolated'],[0,8,16])]:
        for phase in frames:
            data={case:dict(np.load(study/'vortex_structures_v4'/f'{case}_{phase:02d}_sections.npz')) for case in cases}
            for j,span in enumerate(spans):
                key=f'{j}_1_';b,i=data.values();outline=b[f'{j}_actual_outline']
                other=i[f'{j}_actual_outline']
                assert max(cKDTree(outline).query(other)[0].max(),cKDTree(other).query(outline)[0].max())<1e-7
                xx,yy=np.meshgrid(b[key+'x'],b[key+'y']);xc=(xx-outline[:,0].min())/np.ptp(outline[:,0]);gap=yy-top_at(xx,outline)
                region=(xc>=0)&(xc<=.5)&(gap>=.012)&(gap<=.04);mask=region&b[key+'finite']&i[key+'finite']
                row=dict(group=label,phase=phase/16,span=span,coverage=float(mask.sum()/max(region.sum(),1)),cases={})
                for case,d in data.items():
                    q,jump=pressure_jump(d,j);omega=d[key+'omega_span'];a=float(d[key+'spacing'])**2
                    row['cases'][case]=dict(jump=float(np.trapezoid(jump,q)/np.ptp(q)),negative_integral=float(-np.minimum(omega[mask],0).sum()*a/u),
                        rotation_negative=float(-np.minimum(omega[mask&(d[key+'Q']>0)],0).sum()*a/u))
                row['difference']={k:row['cases'][cases[0]][k]-row['cases'][cases[1]][k] for k in row['cases'][cases[0]]};rows.append(row)
    summary=[]
    probes=list(csv.DictReader((study/'section_metrics.csv').open(encoding='utf-8')))
    for span in spans:
        entry=dict(span=span)
        for label,cases in [('moving',['baseline','isolated']),('static',['static_body','static_isolated'])]:
            rs=[r for r in rows if r['span']==span and r['group']==label];t=[r['phase'] for r in rs]
            entry[label]={k:float(np.trapezoid([r['difference'][k] for r in rs],t)) for k in rs[0]['difference']}
            entry[label]['min_coverage']=min(r['coverage'] for r in rs)
            flow=[]
            for case in cases:
                rr=sorted([r for r in probes if r['case']==case and float(r['span'])==span],key=lambda r:float(r['phase']))
                flow.append(float(np.trapezoid([float(r['W_lower_20']) for r in rr],[float(r['phase']) for r in rr])))
            entry[label]['lower_lateral_velocity_difference']=flow[0]-flow[1]
        summary.append(entry)
    out.mkdir(parents=True,exist_ok=True)
    (out/'motion_root_controls.json').write_text(json.dumps(dict(rows=rows,summary=summary,
        definition='Actual upper outline used for each moving/static configuration, matched within each pair. Fixed h=.012–.04L and x/c=0–.5. Difference of moving and static contrasts does not isolate body motion from fin motion.'),indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
