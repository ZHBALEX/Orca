"""Follow the paired layer pressure increment into exact Lamb-vector differences.

Case-midpoint algebra is diagnostic, not a unique causal force partition.
Coordinates, inertial frame, pressure-layer extent and pairing are explicit.
"""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import Delaunay
from run__localPressureBridge import read_profiles,normal_terms
from fn___integratedChain import mean,stats


def curl(gradient):
    return np.stack([gradient[...,2,1]-gradient[...,1,2],gradient[...,0,2]-gradient[...,2,0],gradient[...,1,0]-gradient[...,0,1]],axis=-1)


def quantities(fields,g,u,nu,period,h0,h1):
    b,i,n=[fields[k] for k in ['baseline','isolated','body_only']]
    ds=g['distances'];ix=np.flatnonzero((ds>=h0-1e-9)&(ds<=h1+1e-9));normal=g['normal'][:,:,None,:]
    ub,ui,un=[f['velocity'] for f in [b,i,n]];wb,wi,wn=[curl(f['grad']) for f in [b,i,n]]
    def integrate(vector):return -2*np.trapezoid(np.sum(vector*normal,axis=-1)[:,:,ix],ds[ix],axis=-1)/u**2
    def difference(a,c):return a-c
    q={'lamb':integrate(np.cross(ub,wb)-np.cross(ui,wi)),
       'velocity_change':integrate(np.cross(ub-ui,(wb+wi)/2)),
       'vorticity_change':integrate(np.cross((ub+ui)/2,wb-wi)),
       'body_only_lamb':integrate(np.cross(un,wn))}
    for j,axis in enumerate('xyz'):
        wbj=np.zeros_like(wb);wij=np.zeros_like(wi);wbj[...,j]=wb[...,j];wij[...,j]=wi[...,j]
        q['omega_'+axis]=integrate(np.cross(ub,wbj)-np.cross(ui,wij))
        q['velocity_omega_'+axis]=integrate(np.cross(ub-ui,(wbj+wij)/2))
        q['vorticity_omega_'+axis]=integrate(np.cross((ub+ui)/2,wbj-wij))
    ex=np.zeros_like(normal);ex[...,0]=1
    chord=ex-np.sum(ex*normal,axis=-1,keepdims=True)*normal
    chord/=np.linalg.norm(chord,axis=-1,keepdims=True)
    transverse=np.cross(normal,chord)
    for label,basis in [('chord_tangent',chord),('cross_chord_tangent',transverse),('surface_normal',normal)]:
        wbpart=np.sum(wb*basis,axis=-1,keepdims=True)*basis
        wipart=np.sum(wi*basis,axis=-1,keepdims=True)*basis
        q[label]=integrate(np.cross(ub,wbpart)-np.cross(ui,wipart))
        q[label+'_vorticity_change']=integrate(np.cross((ub+ui)/2,wbpart-wipart))
        q[label+'_velocity_change']=integrate(np.cross(ub-ui,(wbpart+wipart)/2))
    dg=b['grad']-i['grad']
    # Local directions held fixed for each pointwise Cartesian-gradient contraction.
    shear=np.einsum('...i,...ij,...j->...',chord,dg,normal)
    turning=-np.einsum('...i,...ij,...j->...',normal,dg,chord)
    transverse_change=np.sum((wb-wi)*transverse,axis=-1)
    finite=np.isfinite(shear+turning)&np.isfinite(transverse_change)
    np.testing.assert_allclose((shear+turning)[finite],transverse_change[finite],atol=1e-12)
    q['normal_shear_change']=integrate(np.cross((ub+ui)/2,shear[...,None]*transverse))
    q['chord_normal_velocity_change']=integrate(np.cross((ub+ui)/2,turning[...,None]*transverse))
    delta=b['pressure_cp']-i['pressure_cp']
    q.update(surface=b['surface_cp']-i['surface_cp'],outer=delta[:,:,ix[-1]],layer=delta[:,:,ix[0]]-delta[:,:,ix[-1]])
    q['body_only_outer']=n['pressure_cp'][:,:,ix[-1]]
    q['body_only_surface']=n['surface_cp']
    q['body_only_layer']=n['pressure_cp'][:,:,ix[0]]-n['pressure_cp'][:,:,ix[-1]]
    for distance in [.012,.016,.02,.03,.04,.06]:
        index=int(np.flatnonzero(np.isclose(ds,distance))[0]);key=f'{distance:.3f}'
        q['profile_delta_'+key]=delta[:,:,index]
        q['profile_body_'+key]=n['pressure_cp'][:,:,index]
    q['wall']=q['surface']-delta[:,:,ix[0]]
    bt,it=[normal_terms(f,g,ix,u,period,nu,'fourier') for f in [b,i]]
    q.update({key:bt[key]-it[key] for key in bt});q['closure']=q['layer']-sum(q[k] for k in ['lamb','kinetic','unsteady','viscous'])
    for f,v,w in zip([b,i,n],[ub,ui,un],[wb,wi,wn]):
        rebuilt=np.sum(np.cross(v,w)*normal,axis=-1);valid=np.isfinite(rebuilt)&np.isfinite(f['lamb_normal'])
        np.testing.assert_allclose(rebuilt[valid],f['lamb_normal'][valid],atol=1e-12)
    for a,c in [('lamb',q['velocity_change']+q['vorticity_change']),('lamb',sum(q['omega_'+axis] for axis in 'xyz'))]:
        valid=np.isfinite(q[a])&np.isfinite(c);np.testing.assert_allclose(q[a][valid],c[valid],atol=1e-12)
    return q


def paired(q,g,edge,restrict=None):
    c=g['centers'][0];upper=np.flatnonzero(g['upper']);lower=np.flatnonzero(~g['upper'])
    tri=Delaunay(c[lower][:,[0,2]]);query=c[upper][:,[0,2]];simp=tri.find_simplex(query);verts=tri.simplices[simp.clip(0)]
    xy=c[lower][verts][:,:,[0,2]];longest=np.linalg.norm(xy-np.roll(xy,1,axis=1),axis=2).max(1)
    valid=(simp>=0)&(longest<edge);upper=upper[valid];verts=lower[verts[valid]]
    transform=tri.transform[simp[valid]];first=np.einsum('nij,nj->ni',transform[:,:2],query[valid]-transform[:,2]);w=np.column_stack([first,1-first.sum(1)])
    np.testing.assert_allclose(w.sum(1),1,atol=1e-12);assert w.min()>-1e-10
    p={k:np.sum(v[:,verts]*w[None],axis=-1)-v[:,upper] for k,v in q.items()}
    good=np.logical_and.reduce([np.isfinite(v).all(0) for v in p.values()])
    if restrict is not None:good &= np.isin(upper,restrict)
    verts=verts[good];w=w[good];upper=upper[good];p={k:v[:,good] for k,v in p.items()}
    weight=-g['weight'][:,upper];t=g['phase'];signals={k:np.sum(v*weight,axis=1) for k,v in p.items()}
    np.testing.assert_allclose(signals['lamb'],signals['velocity_change']+signals['vorticity_change'],atol=1e-12)
    np.testing.assert_allclose(signals['lamb'],sum(signals['omega_'+axis] for axis in 'xyz'),atol=1e-12)
    np.testing.assert_allclose(signals['lamb'],sum(signals[key] for key in ['chord_tangent','cross_chord_tangent','surface_normal']),atol=1e-12)
    np.testing.assert_allclose(signals['surface_normal'],0,atol=1e-12)
    np.testing.assert_allclose(signals['cross_chord_tangent_vorticity_change'],signals['normal_shear_change']+signals['chord_normal_velocity_change'],atol=1e-12)
    np.testing.assert_allclose(signals['surface'],signals['outer']+signals['layer']+signals['wall'],atol=1e-12)
    rows={k:dict(mean=float(mean(v,t)),quarters=[float(np.trapezoid(v[j:j+5],t[j:j+5])) for j in [0,4,8,12]]) for k,v in signals.items()}
    sides={}
    for key in ['lamb','normal_shear_change','chord_normal_velocity_change','velocity_change','vorticity_change']:
        low=np.sum(np.sum(q[key][:,verts]*w[None],axis=-1)*weight,axis=1)
        up=-np.sum(q[key][:,upper]*weight,axis=1)
        np.testing.assert_allclose(low+up,signals[key],atol=1e-12)
        sides[key]={side:dict(mean=float(mean(value,t)),quarters=[float(np.trapezoid(value[j:j+5],t[j:j+5])) for j in [0,4,8,12]]) for side,value in [('lower_layer',low),('upper_layer',up)]}
    return dict(faces=len(upper),upper_area_coverage=float(g['area'][0,upper].sum()/g['area'][0,g['upper']].sum()),terms=rows,layer_side_contributions=sides),signals,p,upper


def main():
    root=Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    output=Path(__file__).resolve().parents[1]/'reports/2026-10-23-rotation-origin'
    u=.29;nu=1/41379.3;period=1.;edge=.02
    fields,g,manifest=read_profiles(root);records=[];candidates=[];coverage=[]
    for h0,h1 in [(.012,.04),(.008,.03),(.016,.04),(.02,.06)]:
        q=quantities(fields,g,u,nu,period,h0,h1);result,signals,paired_fields,upper=paired(q,g,edge)
        coverage.append(dict(h0=h0,h1=h1,faces=result['faces'],coverage=result['upper_area_coverage'],usable_for_area_comparison=result['upper_area_coverage']>=.9))
        if result['upper_area_coverage']>=.9:candidates.append((h0,h1,q,upper))
    common=candidates[0][3]
    for candidate in candidates[1:]:common=np.intersect1d(common,candidate[3])
    for h0,h1,q,_ in candidates:
        result,signals,paired_fields,upper=paired(q,g,edge,restrict=common)
        result.update(h0=h0,h1=h1);records.append(result)
        if h0==.012:
            main_signals=signals
            shifted={key:value.copy() for key,value in q.items()}
            for key in shifted:
                if key in ['surface','outer'] or key.startswith('profile_delta_'):
                    shifted[key]+=np.sin(2*np.pi*g['phase'])[:,None]+.7
                if key in ['body_only_surface','body_only_outer'] or key.startswith('profile_body_'):
                    shifted[key]+=np.cos(2*np.pi*g['phase'])[:,None]-.4
            _,check_signals,_,_=paired(shifted,g,edge,restrict=common)
            gauge_error=max(float(np.max(abs(check_signals[key]-value))) for key,value in signals.items())
            assert gauge_error<1e-12
            c=g['centers'][0,upper];radius=abs(c[:,2]-6);allradius=abs(g['centers'][0,:,2]-6)
            eta=(radius-allradius.min())/(allradius.max()-allradius.min());bands=[]
            for low,high in [(0,.2),(.2,.4),(.4,.6),(.6,.8),(.8,1.000001)]:
                mask=(eta>=low)&(eta<high);weight=-g['weight'][:,upper[mask]]
                bands.append(dict(lower=low,upper=high,means={k:float(mean(np.sum(v[:,mask]*weight,axis=1),g['phase'])) for k,v in paired_fields.items()}))
            result['span_bands']=bands
    output.mkdir(parents=True,exist_ok=True)
    out=dict(records=records,coverage_audit=coverage,common_faces=len(common),gauge_shift_max_error=gauge_error,provenance=manifest,frame='Saved solver inertial XYZ; x streamwise, y vertical, z lateral. Local chord tangent is global x projected onto the instantaneous surface; cross-chord tangent is normal cross chord tangent, not a reconstructed anatomical span axis.',
        definition='Delta(u cross omega)=Delta u cross midpoint omega + midpoint u cross Delta omega. Midpoints are between cases, not time means.',
        scope='Negative-z fin paired normal layer. Contributions are not a whole-system force or a unique named-vortex force. Velocity/vorticity split is frame-dependent; NoFin Lamb is a control, not a linear causal subtraction.')
    (output/'rotation_origin.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8')
    np.savez_compressed(output/'rotation_signals.npz',phase=g['phase'],**main_signals)
    from fn___plotRotationOrigin import plot_results
    plot_results(output,out)
    print(json.dumps([dict(h0=r['h0'],h1=r['h1'],faces=r['faces'],means={k:r['terms'][k]['mean'] for k in ['lamb','velocity_change','vorticity_change','normal_shear_change','chord_normal_velocity_change','outer','body_only_outer','body_only_surface','body_only_layer']}) for r in records],indent=2),flush=True)


if __name__=='__main__':main()
