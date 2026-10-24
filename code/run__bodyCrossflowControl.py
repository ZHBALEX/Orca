"""Test whether measured lateral-flow changes resemble the body-only field."""
from pathlib import Path
import json
import numpy as np
from run__localPressureBridge import read_profiles


def main():
    root=Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    fields,g,manifest=read_profiles(root);u=.29;t=g['phase'];ix=int(np.flatnonzero(np.isclose(g['distances'],.02))[0])
    radius=abs(g['centers'][0,:,2]-6.);eta=(radius-radius.min())/(radius.max()-radius.min());rows=[]
    for lo,hi in [(0,.2),(.2,.4),(.4,.6),(.6,.8),(.8,1.000001)]:
        for side in ['upper','lower']:
            mask=(eta>=lo)&(eta<hi)&(g['upper'] if side=='upper' else ~g['upper'])
            valid=mask.copy()
            for d in fields.values():valid &= np.isfinite(d['velocity'][:,:,ix,:]).all((0,2))
            weight=g['area'][:,valid];norm=weight.sum(1)
            v={case:np.sum(d['velocity'][:,valid,ix,:]*weight[...,None],axis=1)/norm[:,None]/u for case,d in fields.items()}
            delta=v['baseline']-v['isolated'];body=v['body_only']-np.array([1,0,0]);residual=delta-body
            rows.append(dict(span=[lo,hi],side=side,faces=int(valid.sum()),mean_difference=np.trapezoid(delta,t,axis=0).tolist(),mean_body_disturbance=np.trapezoid(body,t,axis=0).tolist(),
                mean_nonadditive=np.trapezoid(residual,t,axis=0).tolist(),lateral_rms_difference=float(np.sqrt(np.trapezoid(delta[:,2]**2,t))),lateral_rms_nonadditive=float(np.sqrt(np.trapezoid(residual[:,2]**2,t)))))
    out.mkdir(parents=True,exist_ok=True)
    (out/'body_crossflow_control.json').write_text(json.dumps(dict(rows=rows,provenance=manifest,definition='Area-weighted mean Cartesian velocities at matched material normal rays h=.02L on the negative-z fin; body-only sampled at absent-fin fluid positions. Delta velocity compared with NoFin-minus-uniform flow as an additive diagnostic, not a linear governing-equation identity.'),indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
