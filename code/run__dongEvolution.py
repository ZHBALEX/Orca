"""Complete-cycle, matched-window evolution of flow and pressure differences."""
from pathlib import Path
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run__dongComparison import pressure_jump
from fn___plotEvidence import save_hybrid


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-dong-evolution'
    spans=[.2,.3,.4,.5,.7,.9];u=.29
    rows=[]
    for phase in range(17):
        data={case:dict(np.load(study/'vortex_structures_v4'/f'{case}_{phase:02d}_sections.npz')) for case in ['baseline','isolated']}
        for j,span in enumerate(spans):
            key=f'{j}_1_';b,i=data.values()
            xx,yy=np.meshgrid(b[key+'x'],b[key+'y']);xc=(xx-b[key+'leading_x'])/b[key+'chord'];gap=yy-b[key+'upper'][None,:]
            for lo,hi in [(.008,.03),(.012,.04),(.016,.04)]:
                mask=(xc>=0)&(xc<=.5)&(gap>=lo)&(gap<=hi)
                valid=mask&b[key+'finite']&i[key+'finite'];entry=dict(phase=phase/16,span=span,h0=lo,h1=hi,coverage=float(valid.sum()/max(mask.sum(),1)),cases={})
                for case,d in data.items():
                    for name in ['x','y','basis','upper']:
                        np.testing.assert_allclose(d[key+name],b[key+name],atol=1e-10)
                    q,jump=pressure_jump(d,j);omega=d[key+'omega_span'];a=float(d[key+'spacing'])**2
                    entry['cases'][case]=dict(jump=float(np.trapezoid(jump,q)/np.ptp(q)),negative_omega_integral=float(-np.minimum(omega[valid],0).sum()*a/u),
                        rotation_negative_integral=float(-np.minimum(omega[valid&(d[key+'Q']>0)],0).sum()*a/u),
                        span_velocity=float(d[key+'velocity_local'][...,2][valid].mean()/u))
                entry['difference']={k:entry['cases']['baseline'][k]-entry['cases']['isolated'][k] for k in entry['cases']['baseline']}
                rows.append(entry)
    mainrows=[r for r in rows if r['h0']==.012]
    phase=np.arange(17)/16
    keys=['jump','negative_omega_integral','rotation_negative_integral','span_velocity']
    matrices={k:np.array([[next(r for r in mainrows if r['span']==s and r['phase']==t)['difference'][k] for s in spans] for t in phase]) for k in keys}
    summary={}
    for k,m in matrices.items():
        summary[k]=dict(mean_by_span=np.trapezoid(m,phase,axis=0).tolist(),max_phase_by_span=phase[np.argmax(m[:-1],axis=0)].tolist(),min_phase_by_span=phase[np.argmin(m[:-1],axis=0)].tolist())
    # Signed phase relationship is descriptive; frames are not independent samples.
    summary['pressure_vorticity_correlation']=[float(np.corrcoef(matrices['jump'][:-1,j],matrices['negative_omega_integral'][:-1,j])[0,1]) for j in range(6)]
    motion=list(csv.DictReader((study/'section_metrics.csv').open(encoding='utf-8')))
    summary['motion_correlations']=[]
    fits=[]
    for j,span in enumerate(spans):
        selected=sorted([r for r in motion if r['case']=='baseline' and np.isclose(float(r['span']),span)],key=lambda r:float(r['phase']))
        assert len(selected)==17
        summary['motion_correlations'].append(dict(span=span,**{k:float(np.corrcoef(matrices['jump'][:-1,j],np.array([float(r[k]) for r in selected])[:-1])[0,1]) for k in ['body_ay','body_v','pitch_deg']}))
        acceleration=np.array([float(r['body_ay']) for r in selected]);x=-acceleration[:-1];y=matrices['jump'][:-1,j]
        slope,intercept=np.linalg.lstsq(np.column_stack([x,np.ones_like(x)]),y,rcond=None)[0]
        predicted=slope*x+intercept;residual=y-predicted
        fits.append(dict(span=span,slope=float(slope),intercept=float(intercept),R2=float(1-np.sum(residual**2)/np.sum((y-y.mean())**2)),
            residual_vorticity_correlation=float(np.corrcoef(residual,matrices['negative_omega_integral'][:-1,j])[0,1]),
            acceleration=acceleration.tolist(),pressure=matrices['jump'][:,j].tolist()))
    summary['descriptive_acceleration_fits']=fits
    summary['window_sensitivity']=[]
    for lo,hi in [(.008,.03),(.012,.04),(.016,.04)]:
        for span in spans[:3]:
            rs=sorted([r for r in rows if r['h0']==lo and r['span']==span],key=lambda r:r['phase'])
            values=np.array([r['difference']['negative_omega_integral'] for r in rs])
            summary['window_sensitivity'].append(dict(h0=lo,h1=hi,span=span,min_coverage=min(r['coverage'] for r in rs),peak_phase=float(phase[np.argmax(values[:-1])]),mean=float(np.trapezoid(values,phase))))
    out.mkdir(parents=True,exist_ok=True)
    (out/'evolution.json').write_text(json.dumps(dict(rows=rows,summary=summary,scope='Same accepted cycle; local fixed fluid windows, no vortex identity or causality inferred from correlation.'),indent=2),encoding='utf-8')
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(8,5.5),layout='constrained')
    labels=['Pressure-jump increment','Negative-vorticity integral increment / U','Q-positive negative-vorticity increment / U','Spanwise velocity difference / U']
    for ax,k,label in zip(axes.flat,keys,labels):
        for j,color in zip([0,1,2],['#D55E00','#0072B2','#009E73']):
            ax.plot(phase,matrices[k][:,j],'o-',ms=3,color=color,label=f'Span {spans[j]}',zorder=0)
        ax.axhline(0,color='.65',lw=.6,zorder=0);ax.set(xlabel='t/T',ylabel=label);ax.set_rasterization_zorder(1)
        ax.spines[['top','right']].set_visible(False)
    axes[0,0].legend(frameon=False)
    save_hybrid(fig,out/'20_inner_fin_evolution',600)
    fig,axs=plt.subplots(1,3,figsize=(9,3),layout='constrained')
    for ax,f in zip(axs,fits[:3]):
        x=-np.array(f['acceleration'][:-1]);y=np.array(f['pressure'][:-1]);order=np.argsort(x)
        ax.scatter(x,y,c=phase[:-1],cmap='viridis',s=25,zorder=0)
        ax.plot(x[order],(f['slope']*x+f['intercept'])[order],color='black',lw=.7,zorder=0)
        ax.set(xlabel='Negative local vertical acceleration',ylabel='Pressure-jump increment',title=f"Span {f['span']}; R²={f['R2']:.3f}")
        ax.set_rasterization_zorder(1)
    save_hybrid(fig,out/'21_pressure_acceleration',600)
    for svg in out.glob('*.svg'):svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
