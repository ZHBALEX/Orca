"""Phase means at matched chord/clearance coordinates; gauge-free fluid pressure."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import save_hybrid


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-mean-flow-pressure'
    u=.29;section=1;span=.3;phase=np.arange(17)/16
    xc=np.linspace(.05,.9,70);height=np.linspace(.012,.04,24);xx,hh=np.meshgrid(xc,height)
    fields={};sources=[]
    for case in ['baseline','isolated']:
        snapshots=[]
        for k in range(17):
            path=study/'vortex_structures_v4'/f'{case}_{k:02d}_sections.npz';d=np.load(path);prefix=f'{section}_1_'
            x=d[prefix+'x'];y=d[prefix+'y'];physical_x=d[prefix+'leading_x']+xx*d[prefix+'chord']
            top=np.interp(physical_x,x,d[prefix+'upper']);bottom=np.interp(physical_x,x,d[prefix+'lower'])
            upper_query=np.stack([top+hh,physical_x],axis=-1);lower_query=np.stack([bottom-hh,physical_x],axis=-1)
            def sample(key,query):return RegularGridInterpolator((y,x),d[prefix+key],bounds_error=False,fill_value=np.nan)(query)
            v=sample('velocity_local',upper_query);w=sample('omega_span',upper_query)
            jump=2*(sample('pressure',lower_query)-sample('pressure',upper_query))/u**2
            snapshots.append(np.stack([w/u,v[...,2]/u,jump],axis=-1));sources.append(path.name)
        fields[case]=np.stack(snapshots)
    valid=np.isfinite(fields['baseline']).all((0,3)) & np.isfinite(fields['isolated']).all((0,3))
    means={case:np.where(valid[...,None],np.trapezoid(values,phase,axis=0),np.nan) for case,values in fields.items()}
    means['difference']=means['baseline']-means['isolated']
    out.mkdir(parents=True,exist_ok=True)
    metrics=dict(span=span,common_valid_fraction=float(valid.mean()),sampled_coordinates='x/c=.05–.9, equal upper/lower normal-coordinate offsets .012–.04L in each moving local section',
        means={case:[float(np.nanmean(m[...,j])) for j in range(3)] for case,m in means.items()},columns=['upper_span_vorticity_over_U','upper_span_velocity_over_U','lower_minus_upper_fluid_Cp'],
        scope='Unweighted averages over the plotted sampling grid. Not surface pressure or a force integral. Local-coordinate phase means, not stationary Eulerian means.',sources=sources)
    (out/'mean_section.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axs=plt.subplots(3,3,figsize=(9,6.5),layout='constrained')
    for row,label in enumerate(['Upper vorticity / U','Upper span velocity / U','Lower-minus-upper fluid Cp']):
        limit=max(float(np.nanmax(abs(m[...,row]))) for m in means.values())
        for col,(case,m) in enumerate(means.items()):
            ax=axs[row,col];im=ax.pcolormesh(xc,height,m[...,row],cmap='RdBu_r',vmin=-limit,vmax=limit,shading='auto',zorder=0)
            ax.set(xlabel='x/c',ylabel='Clearance h/L',title=case);ax.set_rasterization_zorder(1)
        fig.colorbar(im,ax=axs[row,:],label=label,shrink=.8)
    fig.suptitle('Cycle-mean matched flow and relative pressure at span 0.3')
    save_hybrid(fig,out/'23_mean_flow_pressure',600)
    svg=out/'23_mean_flow_pressure.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':main()
