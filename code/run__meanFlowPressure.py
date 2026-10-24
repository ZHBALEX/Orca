"""Exact cycle-mean product audit on the established paired pressure layer."""
from pathlib import Path
import json
import numpy as np
from run__pairedRotationOrigin import quantities,paired,curl
from run__localPressureBridge import read_profiles


def main():
    root=Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-mean-flow-pressure'
    u,nu,period,edge=.29,1/41379.3,1.,.02
    fields,g,manifest=read_profiles(root)
    common=None
    for h0,h1 in [(.012,.04),(.016,.04),(.02,.06)]:
        upper=paired(quantities(fields,g,u,nu,period,h0,h1),g,edge)[3]
        common=upper if common is None else np.intersect1d(common,upper)
    q=quantities(fields,g,u,nu,period,.012,.04)
    phase=g['phase'];dist=g['distances'];ix=np.flatnonzero((dist>=.012-1e-9)&(dist<=.04+1e-9))
    def avg(v):return np.trapezoid(v,phase,axis=0)
    def integrate(v):return -2*np.trapezoid(np.sum(v*g['normal'][:,:,None,:],axis=-1)[:,:,ix],dist[ix],axis=-1)/u**2
    parts={}
    for case in ['baseline','isolated']:
        vel=fields[case]['velocity'];omega=curl(fields[case]['grad'])
        vm,wm=avg(vel),avg(omega);vp,wp=vel-vm,omega-wm
        parts[case]={'mean_mean':integrate(np.broadcast_to(np.cross(vm,wm),vel.shape)),
            'mean_fluctuation':integrate(np.cross(vm,wp)),
            'fluctuation_mean':integrate(np.cross(vp,wm)),
            'fluctuation_fluctuation':integrate(np.cross(vp,wp))}
        total=integrate(np.cross(vel,omega));reconstructed=sum(parts[case].values());good=np.isfinite(total)&np.isfinite(reconstructed)
        np.testing.assert_allclose(total[good],reconstructed[good],atol=1e-12)
    for key in parts['baseline']:q['lamb_'+key]=parts['baseline'][key]-parts['isolated'][key]
    result,signals,p,upper=paired(q,g,edge,common)
    np.testing.assert_allclose(signals['lamb'],sum(signals['lamb_'+k] for k in parts['baseline']),atol=1e-12)
    radius=abs(g['centers'][0,:,2]-6.);eta=(radius[upper]-radius.min())/(radius.max()-radius.min())
    bands=[]
    for low,high in [(0,.2),(.2,.4),(.4,.6),(.6,.8),(.8,1.000001)]:
        mask=(eta>=low)&(eta<high);weight=-g['weight'][:,upper[mask]]
        bands.append(dict(low=low,high=high,means={k:float(avg((v[:,mask]*weight).sum(1))) for k,v in p.items() if k.startswith('lamb') or k in ['surface','outer','layer','kinetic','unsteady']}))
    out.mkdir(parents=True,exist_ok=True)
    record=dict(result=result,bands=bands,provenance=manifest,
        interpretation='Exact mean/fluctuation algebra at moving material-ray locations in inertial Cartesian components. Mean fields are phase averages at moving sample locations, not Eulerian stationary mean fields. Geometry-dependent projection remains instantaneous; mixed terms need not average to zero. This is a local Lamb pressure-term audit, not whole-system force attribution.')
    (out/'mean_flow_pressure.json').write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in result['terms'].items() if k.startswith('lamb') or k in ['surface','outer','layer','kinetic','unsteady']},indent=2))


if __name__=='__main__':main()
