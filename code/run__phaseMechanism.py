"""Locate pressure enhancement and competing terms within the accepted cycle."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run__pairedRotationOrigin import quantities,paired
from run__localPressureBridge import read_profiles
from fn___plotEvidence import save_hybrid


def main():
    root=Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-23-phase-mechanism'
    u,nu,period,edge=.29,1/41379.3,1.,.02
    fields,g,manifest=read_profiles(root)
    candidates=[quantities(fields,g,u,nu,period,a,b) for a,b in [(.012,.04),(.016,.04),(.02,.06)]]
    common=None
    for q in candidates:
        upper=paired(q,g,edge)[3]
        common=upper if common is None else np.intersect1d(common,upper)
    record,signals,p,upper=paired(candidates[0],g,edge,common)
    radius=abs(g['centers'][0,:,2]-6.)
    eta=(radius[upper]-radius.min())/(radius.max()-radius.min())
    weight=-g['weight'][:,upper];phase=g['phase']
    keys=['surface','outer','wall','layer','lamb','kinetic','unsteady','viscous','closure','normal_shear_change','chord_normal_velocity_change','velocity_change','vorticity_change']
    keys+=['profile_delta_0.060','profile_body_0.060','body_only_outer','body_only_surface']
    maps={k:np.zeros((5,4)) for k in keys}
    for band in range(5):
        mask=(eta>=band/5)&(eta<((band+1)/5 if band<4 else 1.000001))
        for k in keys:
            signal=(p[k][:,mask]*weight[:,mask]).sum(1)
            for quarter,start in enumerate([0,4,8,12]):
                maps[k][band,quarter]=np.trapezoid(signal[start:start+5],phase[start:start+5])
    for k in keys:
        np.testing.assert_allclose(maps[k].sum(),np.trapezoid(signals[k],phase),atol=1e-12)
    np.testing.assert_allclose(maps['surface'],maps['outer']+maps['wall']+maps['layer'],atol=1e-12)
    np.testing.assert_allclose(maps['layer'],sum(maps[k] for k in ['lamb','kinetic','unsteady','viscous','closure']),atol=1e-12)
    order=np.argsort(maps['surface'].ravel())[::-1]
    rows=[]
    for index in order:
        band,quarter=np.unravel_index(index,(5,4))
        rows.append(dict(span_band=int(band),quarter=int(quarter),terms={k:float(v[band,quarter]) for k,v in maps.items()}))
    result=dict(scope='One accepted stable cycle; negative-z common paired fin; identical 1071-face support to v16. Signed integrals over quarter cycles, not conditional quarter means.',
        faces=len(upper),coverage=record['upper_area_coverage'],span_edges=[0,.2,.4,.6,.8,1],phase_edges=[0,.25,.5,.75,1],
        ranked_cells=rows,maps={k:v.tolist() for k,v in maps.items()},means={k:float(v.sum()) for k,v in maps.items()},provenance=manifest)
    out.mkdir(parents=True,exist_ok=True)
    (out/'phase_mechanism.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(10,6),layout='constrained')
    plotted=['surface','outer','lamb','kinetic','unsteady','closure']
    limit=max(abs(maps[k]).max() for k in plotted)*1000
    for ax,k in zip(axes.flat,plotted):
        im=ax.imshow(maps[k]*1000,origin='lower',aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,zorder=0)
        ax.set_rasterization_zorder(1)
        ax.set(title=k.replace('_',' ').capitalize(),xlabel='Phase interval',ylabel='Common-fin span fraction')
        ax.set_xticks(range(4),['0–.25','.25–.50','.50–.75','.75–1'])
        ax.set_yticks(range(5),['0–.2','.2–.4','.4–.6','.6–.8','.8–1'])
    fig.colorbar(im,ax=axes,label='Signed contribution × 1000',shrink=.8)
    save_hybrid(fig,out/'17_phase_pressure_budget',600)
    svg=out/'17_phase_pressure_budget.svg'
    svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(json.dumps(dict(means=result['means'],largest_cells=rows[:4]),indent=2))


if __name__=='__main__':main()
