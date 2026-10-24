"""Audit the new evidence without promoting an unclosed transport budget."""
from pathlib import Path
import json
import numpy as np


def main():
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    read=lambda name:json.loads((out/name).read_text(encoding='utf-8'))
    gates=read('root_gate_flux.json');controls=read('motion_root_controls.json')
    assert len(gates['rows'])==102 and len(controls['rows'])==120
    assert all(abs(c['fluid_gate_fraction']-1)<1e-12 for r in gates['rows'] for c in r['cases'].values())
    assert all(r['coverage']==1 for r in controls['rows'])
    for r in gates['rows']:
        for c in r['cases'].values():np.testing.assert_allclose(c['net_flux'],c['rootward_flux']+c['outward_flux'],atol=1e-14)
    paths=read('root_forward_paths.json')['records'];paired=[]
    for release in [.375,.5]:
        bycase={case:{tuple(r['seed']):r for r in paths if r['case']==case and r['release']==release and r['steps']==1024} for case in ['baseline','isolated']}
        valid=[key for key,r in bycase['baseline'].items() if r['history'] and bycase['isolated'][key]['history']]
        paired.append(dict(release=release,common_valid_seeds=len(valid),approached_upper={case:sum(bycase[case][key]['approached_upper_reference'] for key in valid) for case in bycase}))
    budget=read('upper_vorticity_budget.json')['rows'];transport=[]
    assert len(budget)==90 and min(r['valid_fraction'] for r in budget)==1
    for case in ['baseline','isolated']:
        for span in [.2,.3,.4]:
            rs=[r for r in budget if r['case']==case and r['span']==span]
            ratio=np.sqrt(np.mean([r['rms']['residual']**2 for r in rs]))/np.sqrt(np.mean([r['rms']['time']**2 for r in rs]))
            transport.append(dict(case=case,span=span,residual_rms_over_time_rms=float(ratio)))
    result=dict(gate_count=102,gate_positive_negative_identity='passed',gate_fluid_coverage=1,
        four_case_rows=120,four_case_common_coverage=1,paired_path_control=paired,
        transport_rows=90,transport_coverage=1,transport_residual=transport,
        scientific_verdict='Flow diversion and moving/static pressure contrasts supported. Complete lower-to-upper root leakage suppression not established. Vorticity budget fails quantitative closure; no dominant production pathway or force fraction inferred.',
        gates_remapping='Baseline Cartesian W linearly remapped to isolated-grid gate coordinates; fluid mask nearest-neighbour. Native grids are different.',
        exclusions='Requested gate z=5.925 lacked the defined isolated lower-fin voxel window; omitted, not counted as zero flux. Vorticity time derivatives use only 15 interior native phases at fixed physical positions, not a moving-coordinate time derivative.')
    (out/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
