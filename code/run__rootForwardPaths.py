"""Forward paired particles test the exposed-root route in three dimensions."""
from pathlib import Path
from functools import lru_cache
import sys
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def integrate(seeds,velocity,fluid,start,end,steps):
    dt=1/steps;p=seeds.copy();alive=fluid(start,p)>.999
    histories=[[] for _ in p];reasons=['end' if a else 'invalid_release' for a in alive]
    count=round((end-start)*steps)
    for k in range(count+1):
        t=start+k*dt
        if k%(steps//64)==0 or k==count:
            for j in np.flatnonzero(alive):histories[j].append(dict(t=float(t),xyz=p[j].tolist()))
        if k==count:break
        a=velocity(t,p);p2=p+dt*a/2;b=velocity(t+dt/2,p2);p3=p+dt*b/2;c=velocity(t+dt/2,p3);p4=p+dt*c;d=velocity(t+dt,p4)
        new=p+dt*(a+2*b+2*c+d)/6;valid=np.isfinite(new).all(1)
        for point,time in [(p2,t+dt/2),(p3,t+dt/2),(p4,t+dt),(new,t+dt)]:
            valid &= fluid(time,np.where(np.isfinite(point),point,p))>.999
        for j in np.flatnonzero(alive&~valid):
            reasons[j]='fluid_boundary_or_crop';histories[j].append(dict(t=float(t),xyz=p[j].tolist()))
        p[alive&valid]=new[alive&valid];alive &= valid
    return histories,reasons


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    sys.path.insert(0,str(study/'code_archive/20261023_before_cleanup'))
    from read_volume import nearest_surface
    geom=dict(np.load(study/'loads/isolated.npz'));records=[]
    # Analytic check of the integrator's sign and time interval.
    check,_=integrate(np.zeros((1,3)),lambda t,p:np.ones_like(p),lambda t,p:np.ones(len(p)),.5,1.,512)
    np.testing.assert_allclose(check[0][-1]['xyz'],[.5,.5,.5],atol=1e-12)
    for release in [6,8]:
        ref=dict(np.load(study/'vortex_source_v8/volumes'/f'isolated_{release:02d}.npz'))
        seeds=[];labels=[]
        for x in [10.20,10.23,10.26]:
            for z in [5.905,5.915,5.920]:
                ix=int(np.argmin(abs(ref['x']-x)));iz=int(np.argmin(abs(ref['z']-z)))
                ys=ref['y'][(~ref['fluid'][iz,:,ix].astype(bool))&(ref['y']>5.88)&(ref['y']<6.01)]
                if not len(ys):continue
                for gap in [.006,.012]:
                    seeds.append([float(ref['x'][ix]),float(ys.min()-gap),float(ref['z'][iz])]);labels.append(dict(x=x,z=z,gap=gap))
        for case in ['baseline','isolated']:
            @lru_cache(maxsize=3)
            def frame(j):
                with np.load(study/'vortex_source_v8/volumes'/f'{case}_{j:02d}.npz') as d:
                    axes=(d['z'],d['y'],d['x'])
                    return (RegularGridInterpolator(axes,d['velocity'],bounds_error=False,fill_value=np.nan),RegularGridInterpolator(axes,d['fluid'].astype(float),bounds_error=False,fill_value=0))
            def evaluate(t,p,component):
                tau=np.clip(t,0,1)*16;j=min(int(tau),15);a=tau-j
                return (1-a)*frame(j)[component](p[:,::-1])+a*frame(j+1)[component](p[:,::-1])
            for steps in [512,1024]:
                histories,reasons=integrate(np.array(seeds),lambda t,p:evaluate(t,p,0),lambda t,p:evaluate(t,p,1),release/16,min(1,release/16+.5),steps)
                # Classify proximity against the same isolated-fin reference in both flows.
                indexed={}
                for ident,h in enumerate(histories):
                    for point in h:indexed.setdefault(point['t'],[]).append((ident,point))
                for t,items in indexed.items():
                    j=min(int(t*16),15);a=t*16-j;points=(1-a)*geom['points_history'][j]+a*geom['points_history'][j+1]
                    ids,_,dist=nearest_surface(points,geom['faces'],np.array([p['xyz'] for _,p in items]));vectors=(1-a)*geom['vectors'][j]+a*geom['vectors'][j+1]
                    for n,(_,p) in enumerate(items):
                        p.update(fin_distance=float(dist[n]),nearest_upper=bool(vectors[ids[n],1]<0),reference_face=int(ids[n]))
                for ident,(h,reason) in enumerate(zip(histories,reasons)):
                    approached=any(p['nearest_upper'] and p['fin_distance']<.015 and p['xyz'][0]<10.32 for p in h)
                    records.append(dict(case=case,release=release/16,steps=steps,seed=seeds[ident],label=labels[ident],history=h,reason=reason,approached_upper_reference=approached))
                print(case,release,steps,'done',flush=True)
            frame.cache_clear()
    summary=[]
    for release in [.375,.5]:
        for case in ['baseline','isolated']:
            rs=[r for r in records if r['release']==release and r['case']==case and r['steps']==1024]
            summary.append(dict(release=release,case=case,count=len(rs),valid_release=sum(bool(r['history']) for r in rs),approached_upper=sum(r['approached_upper_reference'] for r in rs)))
    out.mkdir(parents=True,exist_ok=True)
    (out/'root_forward_paths.json').write_text(json.dumps(dict(records=records,summary=summary,scope='Paired lower-side particles. Upper approach means nearest isolated-fin reference face has inward Ay<0, distance<.015L and x<10.32. This proximity is not by itself proof of root bypass; inspect 3D paths. Positions integrate measured 17-phase fields; temporal interpolation is not extra CFD.'),indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
