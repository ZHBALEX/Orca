"""Pressure-to-load accounting with explicit reference and sampling domains."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import Delaunay, cKDTree
from run__pressureJump import paired_chart, jump
from run__localPressureBridge import normal_terms


def mean(v,t):
    return np.trapezoid(v,t,axis=0)/(t[-1]-t[0])


def stats(v,t):
    a=mean(v,t)
    return dict(mean=float(a),rms=float(np.sqrt(mean(v*v,t))),std=float(np.sqrt(mean((v-a)**2,t))),minimum=float(v.min()),maximum=float(v.max()))


def compare(predicted,measured,t):
    error=predicted-measured
    return dict(mean_error=float(mean(error,t)),rms_error=float(np.sqrt(mean(error*error,t))),
        rms_over_measured_fluctuation=float(np.sqrt(mean(error*error,t)/mean((measured-mean(measured,t))**2,t))),
        correlation=float(np.corrcoef(predicted[:-1],measured[:-1])[0,1]))


def system_loads(project,cases,u,area):
    phase=np.arange(481)/480;signals={};rows={};sources=[]
    for name,(folder,start) in cases.items():
        path=project/'Result'/folder/'forces/drag_lift_body_001.dat';data=np.loadtxt(path)
        selected=data[(data[:,0]>=start-1e-5)&(data[:,0]<=start+1+1e-5)]
        assert len(selected)==481,(name,len(selected))
        assert np.max(abs(selected[:,0]-start-phase))<1e-5
        signals[name]={k:np.interp(phase,selected[:,0]-start,selected[:,j])/(u*u*area) for k,j in [('pressure',4),('viscous',5),('total',6)]}
        rows[name]={k:stats(v,phase) for k,v in signals[name].items()}
        sources.append(path)
    signals['interaction']={k:signals['baseline'][k]-signals['isolated'][k]-signals['body_only'][k] for k in signals['baseline']}
    rows['interaction']={k:stats(v,phase) for k,v in signals['interaction'].items()}
    sparse=phase[::30]
    sensitivity={k:dict(dense_mean=float(mean(v,phase)),sparse_mean=float(mean(v[::30],sparse)),
        sparse_mean_error=float(mean(v[::30],sparse)-mean(v,phase))) for k,v in signals['interaction'].items()}
    return dict(case_stats=rows,sampling=sensitivity,
        writer_sum_max_error=float(np.max(abs(signals['interaction']['total']-signals['interaction']['pressure']-signals['interaction']['viscous'])))),phase,signals,sources


def pressure_bridge(data,raw,u,area,edge):
    b,i=data['baseline'],data['isolated'];phase=i['time']-i['time'][0]
    bi,ii,q,vertices,weights,upper=paired_chart(b,i,6.,edge)
    outputs={};all_signals={}
    for name,d,idx in [('baseline',b,bi),('isolated',i,ii)]:
        j=jump(d['marker_pressure'],idx,q,vertices,weights,u)
        a=d['vectors'][:,idx];p=d['marker_pressure'][:,idx];surface_area=np.linalg.norm(a,axis=-1)
        reference=(p[:,~upper]*surface_area[:,~upper]).sum(1)/surface_area[:,~upper].sum(1)
        pair=np.sum(j*(-a[:,q,1]),axis=1)/area
        marker=2*np.sum((p-reference[:,None])*a[:,:,1],axis=1)/(u*u*area)
        raw_p=np.array([r['regions']['common']['pressure_force_y'] for r in raw[name]])
        raw_a=np.array([r['regions']['common']['signed_area_y'] for r in raw[name]])
        native=(raw_p-2*reference*raw_a)/(u*u*area)
        whole=np.array([r['regions']['whole']['pressure_force_y'] for r in raw[name]])/(u*u*area)
        whole_solver=d['solver'][:,4]/(u*u*area)
        # Apply an independent spatially uniform pressure shift to field AND reference.
        shift=.3+np.sin(2*np.pi*phase)*(1 if name=='baseline' else -.7)
        marker_shift=2*np.sum((p+shift[:,None]-(reference+shift)[:,None])*a[:,:,1],axis=1)/(u*u*area)
        native_shift=(raw_p+2*shift*raw_a-2*(reference+shift)*raw_a)/(u*u*area)
        np.testing.assert_allclose(marker_shift,marker,atol=1e-12)
        np.testing.assert_allclose(native_shift,native,atol=1e-12)
        all_signals[name]=dict(pair=pair,marker_reference=marker,native_reference=native,
            geometry_and_unpaired=marker-pair,discretization=native-marker,whole_native=whole,whole_solver=whole_solver)
        outputs[name]=dict(stats={k:stats(v,phase) for k,v in all_signals[name].items()},
            marker_vs_native=compare(marker,native,phase),whole_reconstruction=compare(whole,whole_solver,phase))
    delta={k:all_signals['baseline'][k]-all_signals['isolated'][k] for k in all_signals['baseline']}
    np.testing.assert_allclose(delta['native_reference'],delta['pair']+delta['geometry_and_unpaired']+delta['discretization'],atol=1e-12)
    outputs['delta']=dict(stats={k:stats(v,phase) for k,v in delta.items()},
        pair_vs_native=compare(delta['pair'],delta['native_reference'],phase),marker_vs_native=compare(delta['marker_reference'],delta['native_reference'],phase))
    outputs['reference_definition']='Each case: area-weighted mean pressure on the lower common material surface. Gauge-invariant referenced open-region load, not a uniquely defined physical fin load. Remainder accounts for finite-thickness geometry, unpaired elements and interpolation.'
    outputs['query_faces']=len(q)
    return outputs,phase,delta,all_signals


def profile_budget(fields,g,u,period,nu,h0,h1,method,edge):
    b,i=fields['baseline'],fields['isolated'];ds=g['distances'];ix=np.flatnonzero((ds>=h0-1e-9)&(ds<=h1+1e-9))
    c=g['centers'][0];upper=g['upper'];up=np.flatnonzero(upper);lower=np.flatnonzero(~upper)
    tri=Delaunay(c[lower][:,[0,2]]);query=c[up][:,[0,2]];simplex=tri.find_simplex(query);vertices=tri.simplices[simplex.clip(0)]
    xy=c[lower][vertices][:,:,[0,2]];maximum=np.max(np.linalg.norm(xy-np.roll(xy,1,axis=1),axis=2),axis=1)
    accepted=(simplex>=0)&(maximum<edge);up=up[accepted];v=lower[vertices[accepted]]
    transform=tri.transform[simplex[accepted]];first=np.einsum('nij,nj->ni',transform[:,:2],query[accepted]-transform[:,2]);w=np.column_stack([first,1-first.sum(1)])
    delta=b['pressure_cp']-i['pressure_cp'];q=dict(surface=b['surface_cp']-i['surface_cp'],outer=delta[:,:,ix[-1]],layer=delta[:,:,ix[0]]-delta[:,:,ix[-1]])
    q['wall']=q['surface']-delta[:,:,ix[0]]
    q['lamb']=-2*np.trapezoid(b['lamb_normal'][:,:,ix]-i['lamb_normal'][:,:,ix],ds[ix],axis=-1)/u**2
    bt,it=[normal_terms(p,g,ix,u,period,nu,method) for p in (b,i)]
    q.update({k:bt[k]-it[k] for k in bt})
    q['prediction']=sum(q[k] for k in ['lamb','kinetic','unsteady','viscous']);q['closure']=q['layer']-q['prediction']
    paired={k:np.sum(value[:,v]*w[None],axis=-1)-value[:,up] for k,value in q.items()}
    valid=np.logical_and.reduce([np.isfinite(x).all(0) for x in paired.values()])
    up=up[valid];paired={k:x[:,valid] for k,x in paired.items()}
    np.testing.assert_allclose(paired['surface'],paired['outer']+paired['layer']+paired['wall'],atol=1e-12)
    area=g['area'][:,up];weight=-g['weight'][:,up]
    signals={k:np.sum(value*weight,axis=1) for k,value in paired.items()}
    phase=g['phase'];late=phase>=.75
    result=dict(faces=len(up),common_upper_area_coverage=float(area[0].sum()/g['area'][0,upper].sum()),
        force_stats={k:stats(value,phase) for k,value in signals.items()},
        late_force_means={k:float(mean(value[late],phase[late])) for k,value in signals.items()},
        layer_reconstruction=compare(signals['prediction'],signals['layer'],phase),
        scope='Negative-z fin; jointly valid upper/lower material-chart pairs; paired pressure-jump load proxy, not total fin or whole-system force.')
    return result,signals


def geometry_audit(project,b,post_path):
    import sys
    sys.path.insert(0,str(post_path))
    from tecio_read import read_marker,TecIO
    api=TecIO();folder=project/'Result/NoFin/marker_q';first=read_marker(folder/f'marker{int(b["steps"][0]):09d}.szplt',api=api)
    dist,ids=cKDTree(b['points_history'][0]).query(first['points']);matched=dist<5e-6
    errors=[]
    for j,step in enumerate(b['steps']):
        marker=first if j==0 else read_marker(folder/f'marker{int(step):09d}.szplt',api=api)
        errors.append(float(np.linalg.norm(marker['points'][matched]-b['points_history'][j,ids[matched]],axis=1).max()))
    return dict(matching_tolerance=5e-6,initial_distance_quantiles=np.quantile(dist,[0,.5,.9,1]).tolist(),nofin_vertices=len(dist),matched_vertices=int(matched.sum()),matched_fraction=float(matched.mean()),
        max_matching_motion_error=max(errors),phase_errors=errors,
        limitation='Checks shared vertex trajectories only; unmatched vertices may reflect tessellation or geometry differences; full-surface kinematic equivalence is not established by this check.')


def attribution_audit(study):
    folder=study/'force_partition_v2';checks=[]
    for tag in ['standard','wide','enclosed']:
        for case in ['baseline','isolated']:
            for j in range(17):
                d=json.loads((folder/f'{case}_{j:02d}_{tag}.json').read_text(encoding='utf-8'))
                checks.append(dict(tag=tag,case=case,phase=j/16,error=abs(d['reconstructed_target']-d['native_reference']),potential_residual=d['solver']['relative_residual']))
    prior=json.loads((folder/'attribution_audit_v3.json').read_text(encoding='utf-8'))
    validation=json.loads((folder/'validation_v2.json').read_text(encoding='utf-8'))
    references=prior['provisional_reference']
    selected=next(x for x in references if x['domain']=='wide' and x['excluded_layers']==2) if isinstance(references,list) else references
    d=selected['delta'];target=d['mean_pressure_lift']
    transport=json.loads((study/'vortex_source_v8/resolved_vorticity_transport.json').read_text(encoding='utf-8'))
    windows=[w for r in transport['records'] for w in r.get('windows',[])];valid=[w for w in windows if w.get('gradient_samples_all_valid') and w.get('viscous_samples_all_valid')]
    result=dict(saved_potential_results_checked=len(checks),max_target_identity_error=max(x['error'] for x in checks),
        max_potential_solver_residual=max(x['potential_residual'] for x in checks),
        provisional_delta=d,unassigned_over_target=d['unassigned']/target,
        physical_mean_budget_closed=validation['physical_mean_budget_closed'],
        transport=dict(windows=len(windows),all_terms_valid_windows=len(valid),
            median_abs_residual_over_change=float(np.median([abs(w['residual'])/max(abs(w['omega_x_change']),1e-12) for w in windows])),
            valid_window_median_abs_residual_over_change=float(np.median([abs(w['residual'])/max(abs(w['omega_x_change']),1e-12) for w in valid])),viscosity_provenance=transport['viscosity_provenance']),
        scope='Audit of previously computed potentials and transport windows, not new elliptic solves or new trajectory integration. Old target is a negative-z fin regional load under its original reference convention, not the current closed-system interaction.')
    assert result['max_target_identity_error']<1e-10
    return result
