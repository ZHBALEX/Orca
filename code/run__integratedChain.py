"""Execute the pressure/reference/load/transport evidence chain."""
from pathlib import Path
import json
import hashlib
import numpy as np
from fn___integratedChain import system_loads,pressure_bridge,profile_budget,geometry_audit,attribution_audit,compare
from run__localPressureBridge import read_profiles


def main():
    project=Path(r"D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D")
    study=project/'test/body_fin_mechanism_20261022'
    profiles=Path(r"C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data")
    post=Path(r"C:\Users\ZHB\Documents\GitHub\Picar_Post")
    output=Path(__file__).resolve().parents[1]/'reports/2026-10-23-integrated-chain'
    u=.29;area=.548382;period=1.;nu=1/41379.3;edge=.02
    cases={'baseline':('Sweep35',9.),'isolated':('Sweep35_OnlyFin',3.),'body_only':('NoFin',9.)}
    output.mkdir(parents=True,exist_ok=True)
    raw=json.loads((output/'raw_interface_audit.json').read_text())
    data={}
    for name in ('baseline','isolated'):
        with np.load(study/'loads'/f'{name}.npz') as f:data[name]={k:f[k] for k in f.files}
    system,t,loads,sources=system_loads(project,cases,u,area)
    print('Dense system histories completed',flush=True)
    bridge,phase,delta,regional=pressure_bridge(data,raw,u,area,edge)
    comparison=compare(delta['pair'],loads['interaction']['pressure'][::30],phase)
    fields,g,manifest=read_profiles(profiles)
    budget,signals=profile_budget(fields,g,u,period,nu,.012,.040,'fourier',edge)
    central,_=profile_budget(fields,g,u,period,nu,.012,.040,'central',edge)
    print('Paired normal-momentum budget completed',flush=True)
    geometry=geometry_audit(project,data['baseline'],post)
    attribution=attribution_audit(study)
    sources += [study/'loads'/f'{name}.npz' for name in data]
    sources += sorted((study/'force_partition_v2').glob('*.json'))
    sources += [study/'vortex_source_v8/resolved_vorticity_transport.json',output/'raw_interface_audit.json']
    result=dict(configuration=dict(U=u,Sref=area,T=period,nu=nu,h0=.012,h1=.04,max_pair_edge=edge),
        system=system,pressure_to_load=bridge,pair_vs_system=comparison,paired_normal_budget=budget,
        central_derivative_sensitivity=central,body_geometry=geometry,attribution=attribution,
        provenance=manifest+[dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources])
    (output/'integrated_chain.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    np.savez_compressed(output/'chain_signals.npz',dense_phase=t,phase=phase,system_pressure=loads['interaction']['pressure'],
        system_viscous=loads['interaction']['viscous'],system_total=loads['interaction']['total'],
        **{'regional_'+k:v for k,v in delta.items()},**{'budget_'+k:v for k,v in signals.items()})
    print(json.dumps({k:result[k] for k in ['pair_vs_system','body_geometry','attribution']},indent=2),flush=True)


if __name__=='__main__':main()
