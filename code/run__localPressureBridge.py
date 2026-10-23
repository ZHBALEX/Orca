"""Exploratory material-patch pressure test near the previously selected RV+ core."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import style_axes, save_hybrid
from fn___auditEvidence import time_mean


def read_profiles(root):
    fields={}; manifest=[]
    for case in ('baseline','isolated','body_only'):
        frames=[]
        for j in range(17):
            path=root/f'{case}_{j:02d}_profiles.npz'
            with np.load(path) as f:
                frames.append({k:f[k] for k in ('surface_cp','pressure_cp','lamb_normal','velocity','grad','kinetic_gradient_normal','laplacian_normal')})
            manifest.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        fields[case]={k:np.stack([f[k] for f in frames]) for k in frames[0]}
    with np.load(root/'geometry.npz') as f:
        geometry={k:f[k] for k in ('centers','upper','area','weight','distances','phase','normal','ray_velocity')}
    manifest.append(dict(file='geometry.npz',sha256=hashlib.sha256((root/'geometry.npz').read_bytes()).hexdigest()))
    return fields,geometry,manifest


def normal_terms(profile,g,ix,u,period,nu,method):
    first=profile['velocity'][:16]
    if method=='central':
        derivative=(np.roll(first,-1,axis=0)-np.roll(first,1,axis=0))*8/period
    else:
        frequency=np.fft.fftfreq(16,d=period/16);frequency[8]=0
        derivative=np.fft.ifft(np.fft.fft(first,axis=0)*(2j*np.pi*frequency)[:,None,None,None],axis=0).real
    derivative=np.concatenate([derivative,derivative[:1]])
    eulerian=derivative-np.einsum('...ij,...j->...i',profile['grad'],g['ray_velocity'])
    normal=g['normal'][:,:,None,:]
    integrands=dict(unsteady=np.sum(eulerian*normal,axis=-1),kinetic=profile['kinetic_gradient_normal'],viscous=-nu*profile['laplacian_normal'])
    return {k:2*np.trapezoid(v[:,:,ix],g['distances'][ix],axis=-1)/u**2 for k,v in integrands.items()}


def measure(fields,g,core,anchor_frame,radii,h0,h1,u,period,nu,method):
    b,i,no=[fields[k] for k in ('baseline','isolated','body_only')]
    ds=g['distances']; ix=np.flatnonzero((ds>=h0-1e-9)&(ds<=h1+1e-9))
    delta=b['surface_cp']-i['surface_cp']
    incident=no['surface_cp']
    residual=delta-incident
    lamb=-2*np.trapezoid(b['lamb_normal'][:,:,ix]-i['lamb_normal'][:,:,ix],ds[ix],axis=-1)/u**2
    lamb_residual=-2*np.trapezoid(b['lamb_normal'][:,:,ix]-i['lamb_normal'][:,:,ix]-no['lamb_normal'][:,:,ix],ds[ix],axis=-1)/u**2
    residual_profile=b['pressure_cp']-i['pressure_cp']-no['pressure_cp']
    outer=residual_profile[:,:,ix[-1]]
    layer=residual_profile[:,:,ix[0]]-outer
    wall=residual-residual_profile[:,:,ix[0]]
    quantities=dict(delta=delta,incident=incident,residual=residual,lamb=lamb,
                    lamb_residual=lamb_residual,outer_residual=outer,layer_residual=layer,wall_residual=wall)
    bt,it,nt=[normal_terms(p,g,ix,u,period,nu,method) for p in (b,i,no)]
    for key in bt:quantities[key+'_residual']=bt[key]-it[key]-nt[key]
    quantities['predicted_residual']=sum(quantities[k+'_residual'] for k in ('lamb','unsteady','kinetic','viscous'))
    quantities['closure_residual']=layer-quantities['predicted_residual']
    valid=g['upper'].copy()
    for v in quantities.values(): valid &= np.isfinite(v).all(0)
    np.testing.assert_allclose(residual[:,valid],(outer+layer+wall)[:,valid],atol=1e-12)
    center_distance=np.linalg.norm(g['centers'][anchor_frame]-core,axis=1)
    nearest=int(np.argmin(np.where(g['upper'],center_distance,np.inf)))
    anchor=g['centers'][anchor_frame,nearest]
    surface_distance=np.linalg.norm(g['centers'][anchor_frame]-anchor,axis=1)
    masks={'all_upper':valid}
    masks.update({f'radius_{radius:.3f}':valid&(surface_distance<=radius) for radius in radii})
    rows=[];signals={}
    for label,mask in masks.items():
        if not np.any(mask): continue
        area=g['area'][:,mask]; weight=g['weight'][:,mask]
        s={}
        for key,v in quantities.items():
            s[key+'_cp']=np.sum(v[:,mask]*area,axis=1)/area.sum(1)
            s[key+'_cl']=np.sum(v[:,mask]*weight,axis=1)
        for kind in ('cp','cl'):
            np.testing.assert_allclose(s['residual_'+kind],sum(s[k+'_residual_'+kind] for k in ('outer','layer','wall')),atol=1e-12)
        phase=g['phase']; late=phase>=.75
        rows.append(dict(patch=label,faces=int(mask.sum()),area_fraction=float(area[anchor_frame].sum()/g['area'][anchor_frame,valid].sum()),
                         cycle_mean={k:time_mean(v,phase) for k,v in s.items()},
                         late_window_mean={k:time_mean(v[late],phase[late]) for k,v in s.items()},
                         late_window_integral={k:float(np.trapezoid(v[late],phase[late])) for k,v in s.items()}))
        signals[label]=s
    total=rows[0]
    for row in rows:
        row['late_residual_lift_share']=row['late_window_integral']['residual_cl']/total['late_window_integral']['residual_cl']
        row['cycle_residual_lift_share']=row['cycle_mean']['residual_cl']/total['cycle_mean']['residual_cl']
    assert all(rows[j]['faces']<rows[j+1]['faces'] for j in range(1,len(rows)-1))
    return dict(rows=rows,anchor_surface_face=nearest,anchor_surface_xyz=anchor.tolist(),
                core_to_upper_face_distance=float(center_distance[nearest]),valid_upper_faces=int(valid.sum()),
                valid_upper_area_fraction=float(g['area'][anchor_frame,valid].sum()/g['area'][anchor_frame,g['upper']].sum()),
                checks='Surface=outer+layer+wall verified per face and in each patch pressure/force signal; nested patch counts verified.',
                interpretation='Exploratory same-material spatial association, NOT an RV+ force fraction. Surface distance is Euclidean at phase 0.75; patches remain fixed by material index. Body-only pressure is an incident-field control, not a unique physical partition.'),signals,masks


def plot_bridge(g,result,signals,masks,origin,length,colors):
    fig,axs=plt.subplots(1,3,figsize=(10,3.4))
    fig.subplots_adjust(left=.07,right=.99,bottom=.22,top=.86,wspace=.47)
    for ax,p in zip(axs,'abc'):style_axes(ax,p)
    frame=12; xyz=(g['centers'][frame]-origin)/length
    axs[0].scatter(xyz[g['upper'],0],xyz[g['upper'],2],s=2,c='#cccccc',zorder=0)
    names=[r['patch'] for r in result['rows'][1:]]
    for name,c in reversed(list(zip(names,colors))):
        m=masks[name];axs[0].scatter(xyz[m,0],xyz[m,2],s=5,c=c,zorder=0)
    axs[0].set(xlabel=r'$x^*$',ylabel=r'$z^*$',title='Material patches',aspect='equal',ylim=(-.23,-.075))
    for name,c in zip(names,colors):
        axs[1].plot(g['phase'],signals[name]['residual_cp'],'o-',ms=2,lw=1,color=c,zorder=0,label=name.replace('radius_','r/L='))
    axs[1].plot(g['phase'],signals['all_upper']['residual_cp'],'k--',lw=1,zorder=0,label='All valid upper')
    axs[1].set(xlabel=r'$t/T$',ylabel=r'$\Delta C_p-C_{p,\mathrm{NoFin}}$',title='Incident-removed surface pressure')
    axs[1].legend(frameon=False,fontsize=7,loc='best')
    rows=result['rows'][1:];x=np.arange(len(rows))
    axs[2].bar(x-.18,[100*r['area_fraction'] for r in rows],width=.36,color='#92979b',label='Area',zorder=0)
    axs[2].bar(x+.18,[100*r['late_residual_lift_share'] for r in rows],width=.36,color='#d71920',label='Residual lift',zorder=0)
    axs[2].set_xticks(x,[n.replace('radius_','') for n in names]);axs[2].set(xlabel=r'$r/L$',ylabel='Share (%)',title=r'Late window: $0.75\leq t/T\leq1$')
    axs[2].legend(frameon=False,fontsize=8)
    return fig


def plot_local_terms(result):
    fig,axs=plt.subplots(1,2,figsize=(7.5,3.5))
    fig.subplots_adjust(left=.11,right=.98,bottom=.30,top=.83,wspace=.38)
    keys=['lamb_residual_cp','kinetic_residual_cp','unsteady_residual_cp','viscous_residual_cp','predicted_residual_cp','layer_residual_cp']
    for ax,p,row in zip(axs,'ab',[result['rows'][1],result['rows'][-1]]):
        style_axes(ax,p)
        ax.bar(range(6),[row['late_window_mean'][k] for k in keys],color=['#d71920','#087bc1','#a11c68','#92979b','#666666','#087bc1'],zorder=0)
        ax.set_xticks(range(6),['Lamb','Kinetic','Unsteady','Viscous','Rebuilt','Measured'],rotation=40,ha='right')
        ax.set(ylabel='Incident-removed layer pressure',title=row['patch'].replace('radius_',r'$r/L=$'),ylim=(-.22,.09))
    return fig


def main():
    profile_root=Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    output=Path(__file__).resolve().parents[1]/'reports'/'2026-10-23-evidence-audit'
    core=np.array([10.273124694824219,5.939525127410889,5.912974834442139])
    anchor_frame=12
    radii=[.015,.025,.040]
    h0,h1=.012,.040
    u,length=.29,1.
    period,nu=1.,1/41379.3
    method='fourier'
    origin=np.array([10.,6.,6.])
    colors=['#d71920','#087bc1','#a11c68']
    dpi=600
    matplotlib.rcParams.update({'font.family':'Times New Roman','font.size':10,'mathtext.fontset':'stix','svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7})
    output.mkdir(parents=True,exist_ok=True)
    fields,g,manifest=read_profiles(profile_root)
    result,signals,masks=measure(fields,g,core,anchor_frame,radii,h0,h1,u,period,nu,method)
    central,_,_=measure(fields,g,core,anchor_frame,radii,h0,h1,u,period,nu,'central')
    result['central_method_rows']=central['rows']
    result['source_manifest']=manifest
    result['settings']=dict(core_xyz=core.tolist(),anchor_frame=anchor_frame,radii=radii,h0=h0,h1=h1,U=u,L=length,period=period,nu=nu,method=method)
    result['figure_structure']=save_hybrid(plot_bridge(g,result,signals,masks,origin,length,colors),output/'03_local_pressure_bridge',dpi)
    result['local_terms_figure_structure']=save_hybrid(plot_local_terms(result),output/'04_local_momentum_terms',dpi)
    (output/'local_pressure_bridge.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='source_manifest'},indent=2))


if __name__=='__main__':main()
