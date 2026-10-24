"""Discriminate geometric, static-body, incidence and motion-related explanations.

Exploratory analysis of existing four-case simulations, not a preregistered test.
The pressure estimand is the same material-paired common-fin load proxy in all cases.
"""
from pathlib import Path
import csv
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run__pressureJump import paired_chart,jump,integrate_time as avg
from fn___plotEvidence import save_hybrid,style_axes


def controls(data,u,area):
    t=np.arange(17)/16
    bi,ii,q,v,w,upper=paired_chart(data['baseline'],data['isolated'],6.,.02)
    ids={'baseline':bi,'isolated':ii,'static_body':bi,'static_isolated':ii}
    checks={};signals={};jumps={};weights={};case_stats={}
    for dyn,sta in [('baseline','static_body'),('isolated','static_isolated')]:
        np.testing.assert_array_equal(data[dyn]['faces'],data[sta]['faces'])
        np.testing.assert_array_equal(data[dyn]['common_mask'],data[sta]['common_mask'])
        checks[dyn+'_static_mean_shape_error']=float(np.max(abs(avg(data[dyn]['centers'],t)-data[sta]['centers'][0])))
        assert checks[dyn+'_static_mean_shape_error']<1e-6
    np.testing.assert_allclose(data['static_body']['centers'][:,bi],data['static_isolated']['centers'][:,ii],atol=1e-7)
    np.testing.assert_allclose(data['baseline']['vectors'][:,bi],data['isolated']['vectors'][:,ii],atol=1e-10)
    for name,d in data.items():
        idx=ids[name];j=jump(d['marker_pressure'],idx,q,v,w,u);weight=-d['vectors'][:,idx[q],1]/area
        shifted=jump(d['marker_pressure']+(.2+np.sin(2*np.pi*t))[:,None],idx,q,v,w,u)
        np.testing.assert_allclose(j,shifted,atol=1e-12)
        jumps[name]=j;weights[name]=weight;signals[name]=np.sum(j*weight,axis=1)
        case_stats[name]=dict(mean=float(avg(signals[name],t)),std=float(np.sqrt(avg((signals[name]-avg(signals[name],t))**2,t))),
            max_surface_motion=float(np.max(np.linalg.norm(d['centers']-d['centers'][0],axis=-1))))
    dynamic=jumps['baseline']-jumps['isolated'];static=jumps['static_body']-jumps['static_isolated']
    signals['dynamic_contrast']=signals['baseline']-signals['isolated']
    signals['static_contrast']=signals['static_body']-signals['static_isolated']
    signals['difference_in_differences']=signals['dynamic_contrast']-signals['static_contrast']
    dm=float(avg(signals['dynamic_contrast'],t));sm=float(avg(signals['static_contrast'],t))
    weight=weights['isolated'];mean_pressure=float(np.sum(avg(dynamic,t)*avg(weight,t)))
    covariance=dm-mean_pressure
    # Alternate chart built on the static geometry tests the pairing convention.
    sb,si,sq,sv,sw,_=paired_chart(data['static_body'],data['static_isolated'],6.,.02)
    sj=jump(data['static_body']['marker_pressure'],sb,sq,sv,sw,u)-jump(data['static_isolated']['marker_pressure'],si,sq,sv,sw,u)
    static_chart=float(avg((sj*(-data['static_isolated']['vectors'][:,si[sq],1])/area).sum(1),t))
    c=data['isolated']['centers'][0,ii[q]];distance=abs(c[:,2]-6.)
    common_centers=data['isolated']['centers'][0,ii];radius=abs(common_centers[:,2]-6.)
    eta=(distance-radius.min())/(radius.max()-radius.min());edges=np.array([-np.inf,.2,.4,.6,.8,np.inf]);bands=[]
    for a,b in zip(edges[:-1],edges[1:]):
        mask=(eta>=a)&(eta<b);sig=np.sum(dynamic[:,mask]*weight[:,mask],axis=1)
        bands.append(dict(mean=float(avg(sig,t)),area=float(np.linalg.norm(data['isolated']['vectors'][0,ii[q[mask]]],axis=-1).sum()),
            quarter_integrals=[float(np.trapezoid(sig[k:k+5],t[k:k+5])) for k in [0,4,8,12]]))
    np.testing.assert_allclose(sum(x['mean'] for x in bands),dm,atol=1e-12)
    result=dict(case_stats=case_stats,checks=checks,paired_faces=len(q),dynamic_mean=dm,static_mean=sm,
        static_fraction_of_dynamic=sm/dm,motion_associated_excess=dm-sm,motion_associated_fraction=(dm-sm)/dm,
        alternate_static_chart_mean=static_chart,alternate_static_chart_faces=len(sq),
        mean_pressure_on_mean_projection=mean_pressure,pressure_projection_covariance=covariance,
        projection_covariance_fraction=covariance/dm,bands=bands,span_definition='Five equal bands of the common material-surface span, root to tip.',
        scope='Paired common-fin pressure load proxy; not full system lift. Motion contrast freezes the entire geometry at its mean pose, so it does not isolate body motion from fin motion.')
    return result,t,signals


def inflow(path):
    rows=list(csv.DictReader(path.open(encoding='utf-8')));results=[];profile_audit=[]
    for span in sorted({float(r['span']) for r in rows}):
        cases={}
        for name in ['baseline','isolated','static_body','static_isolated']:
            selected=sorted([r for r in rows if r['case']==name and float(r['span'])==span],key=lambda r:float(r['phase']))
            t=np.array([float(r['phase']) for r in selected]);assert len(t)==(3 if name.startswith('static') else 17)
            keys=['Cp_lower','Cp_upper','body_ay','body_v','pitch_deg']+[f'{k}_{distance}' for distance in [15,30,50] for k in ['alpha','speed2','qs_proxy']]
            values={k:np.array([float(r[k]) for r in selected]) for k in keys}
            assert all(np.isfinite(x).all() for x in values.values())
            section_index=sorted({float(r['span']) for r in rows}).index(span)
            for r in selected:
                cache=path.parent/'sections'/f'{name}_{round(float(r["phase"])*16):02d}.npz'
                with np.load(cache) as f:
                    query=f[f'{section_index}_query'];profiles=f[f'{section_index}_profiles']
                    measured=np.trapezoid(profiles[1]-profiles[0],query)/(query[-1]-query[0])
                error=abs(measured-(float(r['Cp_lower'])-float(r['Cp_upper'])))
                assert error<1e-12
                profile_audit.append(float(error))
            cases[name]={k:float(avg(x,t)) for k,x in values.items()}
            cases[name]['sample_count']=len(t)
            cases[name]['pressure_jump']=float(avg(values['Cp_lower']-values['Cp_upper'],t))
        b,i=cases['baseline'],cases['isolated'];results.append(dict(span=span,cases=cases,
            delta_jump=b['pressure_jump']-i['pressure_jump'],probe_comparison={str(d):dict(delta_alpha=b[f'alpha_{d}']-i[f'alpha_{d}'],
            relative_speed_squared_ratio=b[f'speed2_{d}']/i[f'speed2_{d}'],delta_qs_proxy=b[f'qs_proxy_{d}']-i[f'qs_proxy_{d}']) for d in [15,30,50]}))
    return results,dict(section_phase_profiles_checked=len(profile_audit),max_jump_reintegration_error=max(profile_audit))


def plots(output,result,t,signals):
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    blue='#087bc1';red='#d71920';grey='#92979b'
    fig,axs=plt.subplots(1,3,figsize=(10.8,3.8));fig.subplots_adjust(left=.07,right=.94,bottom=.26,top=.84,wspace=.45)
    for ax,label in zip(axs,'abc'):style_axes(ax,label)
    c=result['controls'];names=['static_isolated','static_body','isolated','baseline']
    axs[0].bar(range(4),[1000*c['case_stats'][k]['mean'] for k in names],color=[grey,blue,grey,red],zorder=0)
    axs[0].set_xticks(range(4),['Static I','Static B','Moving I','Moving B'],rotation=35,ha='right')
    axs[0].set(ylabel=r'$10^3\overline{C}_{L,p}^{\rm pair}$',title='Four existing simulations')
    axs[1].bar(range(3),[1000*c[k] for k in ['static_mean','motion_associated_excess','dynamic_mean']],color=[blue,red,grey],zorder=0)
    axs[1].set_xticks(range(3),['Static B-I','Motion excess','Moving B-I'],rotation=35,ha='right')
    axs[1].set(ylabel=r'$10^3\overline{\Delta C}_{L,p}^{\rm pair}$',title='Effect of freezing motion')
    matrix=np.array([b['quarter_integrals'] for b in c['bands']])*1000
    axs[2].lines[0].remove()  # The heat-map row index is not a physical zero axis.
    im=axs[2].imshow(matrix,aspect='auto',cmap='RdBu_r',vmin=-abs(matrix).max(),vmax=abs(matrix).max(),zorder=0)
    for row in range(5):
        for col in range(4):axs[2].text(col,row,f'{matrix[row,col]:.2f}',ha='center',va='center',fontsize=7,color='white' if abs(matrix[row,col])>.7*abs(matrix).max() else 'black')
    axs[2].set_xticks(range(4),['0–.25','.25–.5','.5–.75','.75–1'],rotation=35,ha='right');axs[2].set_yticks(range(5),['0–.2','.2–.4','.4–.6','.6–.8','.8–1'])
    axs[2].set(xlabel=r'$t/T$',ylabel='Common span (root to tip)',title='Signed quarter contributions')
    fig.colorbar(im,ax=axs[2],fraction=.045,pad=.03)
    save_hybrid(fig,output/'12_motion_controls',600)
    fig,axs=plt.subplots(1,3,figsize=(10.8,3.5));fig.subplots_adjust(left=.075,right=.99,bottom=.20,top=.84,wspace=.45)
    for ax,label in zip(axs,'abc'):style_axes(ax,label)
    rows=result['inflow'];span=[r['span'] for r in rows]
    for distance,color in [('15',red),('30',blue),('50',grey)]:
        axs[0].plot(span,[r['probe_comparison'][distance]['delta_alpha'] for r in rows],'o-',color=color,ms=3,lw=1,zorder=0,label=distance+'% chord upstream')
        axs[1].plot(span,[100*(r['probe_comparison'][distance]['relative_speed_squared_ratio']-1) for r in rows],'o-',color=color,ms=3,lw=1,zorder=0)
    axs[0].set(xlabel='Span fraction',ylabel=r'$\overline{\alpha}_B-\overline{\alpha}_I$ (deg)',title='Incidence does not rise inboard')
    axs[0].legend(frameon=False,fontsize=7)
    axs[1].set(xlabel='Span fraction',ylabel='Change in mean relative speed squared (%)',title='Relative inflow changes')
    for pair,label,color in [(('baseline','isolated'),'Moving',red),(('static_body','static_isolated'),'Static',blue)]:
        axs[2].plot(span,[r['cases'][pair[0]]['pressure_jump']-r['cases'][pair[1]]['pressure_jump'] for r in rows],'o-',color=color,ms=3,lw=1,zorder=0,label=label)
    axs[2].set(xlabel='Span fraction',ylabel=r'$\overline{\Delta J}_{\rm section}$',title='Mean pressure-jump response')
    axs[2].legend(frameon=False,fontsize=8)
    save_hybrid(fig,output/'13_inflow_and_pressure',600)
    for path in output.glob('*.svg'):path.write_text('\n'.join(x.rstrip() for x in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')


def source_terms(folder):
    rows=[];paths=[];t=np.arange(17)/16
    for tag in ['standard','wide','enclosed']:
        records={}
        for case in ['baseline','isolated']:
            records[case]=[]
            for j in range(17):
                path=folder/f'{case}_{j:02d}_{tag}.json';paths.append(path)
                records[case].append(json.loads(path.read_text(encoding='utf-8')))
        mean=lambda key:float(avg(np.array([a[key]-b[key] for a,b in zip(records['baseline'],records['isolated'])]),t))
        row=dict(domain=tag,acceleration=mean('acceleration'),outer=mean('outer'),original_regional_target=mean('native_reference'),Q_sensitivity=[])
        if tag!='enclosed':
            for erosion in [1,2,3]:
                def qm(key):
                    values={case:[next(q for q in r['Q'] if q['erosion']==erosion)[key] for r in rs] for case,rs in records.items()}
                    return float(avg(np.array(values['baseline'])-values['isolated'],t))
                q={k:qm(k) for k in ['force','rotation','strain','unassigned']};q['erosion']=erosion
                np.testing.assert_allclose(q['force'],q['rotation']+q['strain'],atol=1e-12)
                np.testing.assert_allclose(row['original_regional_target'],row['acceleration']+row['outer']+q['force']+q['unassigned'],atol=1e-10)
                row['Q_sensitivity'].append(q)
        rows.append(row)
    return dict(rows=rows,scope='Reintegrated existing auxiliary-potential results. Original one-fin target/reference; do not combine its numerical fractions with paired-fin controls. Acceleration and Q source terms are pressure-gauge independent; outer and unassigned bookkeeping retains the original convention.'),paths


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    output=Path(__file__).resolve().parents[1]/'reports/2026-10-23-lift-cause-controls'
    u=.29;area=.548382
    data={};paths=[]
    for name in ['baseline','isolated','static_body','static_isolated']:
        path=study/'loads'/f'{name}.npz';paths.append(path)
        with np.load(path) as f:data[name]={k:f[k] for k in f.files}
    result,t,signals=controls(data,u,area)
    path=study/'section_metrics.csv';paths.append(path)
    output.mkdir(parents=True,exist_ok=True)
    inflow_rows,profile_audit=inflow(path)
    source_result,source_paths=source_terms(study/'force_partition_v2');paths+=source_paths
    paths += [p for p in sorted((study/'sections').glob('*.npz')) if p.stem.rsplit('_',1)[0] in data]
    metrics=dict(controls=result,inflow=inflow_rows,profile_audit=profile_audit,source_terms=source_result,provenance=[dict(name=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
        declared_rivals=['Pressure-zero artifact','Different common-fin motion or mean shape','Pressure/projection covariance','Static-body effect alone','Larger mean incidence alone','Motion-dependent body-fin flow response'],
        interpretation='Post hoc discriminating controls, not unique vortex-force attribution. Upstream probes include self-induced flow and are not independent imposed inflow.')
    (output/'cause_controls.json').write_text(json.dumps(metrics,indent=2,allow_nan=False),encoding='utf-8')
    np.savez_compressed(output/'cause_signals.npz',phase=t,**signals)
    plots(output,metrics,t,signals)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
