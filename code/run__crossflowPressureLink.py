"""Check whether lateral velocity change directly supplies the local pressure term.

Exact Cartesian midpoint algebra; not a counterfactual CFD solution.
"""
from pathlib import Path
import json
import numpy as np
from run__pairedRotationOrigin import quantities,paired,curl
from run__localPressureBridge import read_profiles


def main():
    root=Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    fields,g,manifest=read_profiles(root);u=.29;common=None
    for h0,h1 in [(.012,.04),(.016,.04),(.02,.06)]:
        upper=paired(quantities(fields,g,u,1/41379.3,1,h0,h1),g,.02)[3]
        common=upper if common is None else np.intersect1d(common,upper)
    q=quantities(fields,g,u,1/41379.3,1,.012,.04);ds=g['distances'];ix=np.flatnonzero((ds>=.012-1e-9)&(ds<=.04+1e-9))
    ub,ui=[fields[c]['velocity'] for c in ['baseline','isolated']]
    wb,wi=[curl(fields[c]['grad']) for c in ['baseline','isolated']]
    du=ub-ui;wm=(wb+wi)/2;um=(ub+ui)/2;dw=wb-wi
    def integral(v):return -2*np.trapezoid(np.sum(v*g['normal'][:,:,None,:],axis=-1)[:,:,ix],ds[ix],axis=-1)/u**2
    for j,axis in enumerate('xyz'):
        delta=np.zeros_like(du);delta[...,j]=du[...,j]
        q['direct_velocity_'+axis]=integral(np.cross(delta,wm))
        midpoint=np.zeros_like(um);midpoint[...,j]=um[...,j]
        q['changed_vorticity_with_velocity_'+axis]=integral(np.cross(midpoint,dw))
    result,signals,p,upper=paired(q,g,.02,common)
    np.testing.assert_allclose(signals['velocity_change'],sum(signals['direct_velocity_'+a] for a in 'xyz'),atol=1e-12)
    np.testing.assert_allclose(signals['vorticity_change'],sum(signals['changed_vorticity_with_velocity_'+a] for a in 'xyz'),atol=1e-12)
    (out/'crossflow_pressure_link.json').write_text(json.dumps(dict(result=result,provenance=manifest,
        scope='Exact local paired-layer algebra in solver inertial XYZ. Delta(u cross omega)=Delta u cross midpoint omega + midpoint u cross Delta omega. Separates direct velocity effects from changed-vorticity effects, but does not identify the cause of Delta omega. Not a modified-flow simulation or whole-force partition.'),indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v['mean'] for k,v in result['terms'].items() if 'velocity_' in k or k in ['lamb','vorticity_change']},indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from fn___plotEvidence import save_hybrid
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axs=plt.subplots(1,2,figsize=(8,3.4),layout='constrained')
    for ax,prefix,title in zip(axs,['direct_velocity_','changed_vorticity_with_velocity_'],['Velocity change × midpoint vorticity','Midpoint velocity × vorticity change']):
        ax.bar(range(3),[1000*result['terms'][prefix+a]['mean'] for a in 'xyz'],color=['#0072B2','#009E73','#D55E00'],zorder=0)
        ax.set_xticks(range(3),['Streamwise U','Vertical V','Lateral W'])
        ax.set(title=title,ylabel='Paired-layer pressure term × 1000',ylim=(-.5,2.6))
        ax.axhline(0,color='.5',lw=.6,zorder=0);ax.set_rasterization_zorder(1)
        ax.spines[['top','right']].set_visible(False)
    save_hybrid(fig,out/'27_velocity_vorticity_pressure_link',600)
    svg=out/'27_velocity_vorticity_pressure_link.svg'
    svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')


if __name__=='__main__':main()
