"""Gauge-invariant, material-chart upper/lower pressure comparison."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree, Delaunay
from fn___spatialEvidence import match_area
from fn___plotEvidence import save_hybrid, style_axes


def read_loads(root):
    data = {}; manifest = []
    for name in ('baseline', 'isolated'):
        path = root/f'{name}.npz'
        with np.load(path) as f:
            data[name] = {k: f[k] for k in ('centers', 'vectors', 'marker_pressure', 'common_mask', 'pectoral_mask', 'time')}
        manifest.append(dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return data, manifest


def paired_chart(b, i, midplane, max_edge):
    bi = np.flatnonzero(b['common_mask']); ii = np.flatnonzero(i['common_mask'])
    distance, index = cKDTree(b['centers'][0, bi]).query(i['centers'][0, ii])
    assert np.max(distance) < 1e-8 and len(np.unique(index)) == len(index)
    bi = bi[index]
    np.testing.assert_allclose(b['centers'][:, bi], i['centers'][:, ii], atol=1e-7)
    c = i['centers'][0, ii]; v = i['vectors'][0, ii]
    upper = v[:, 1] < 0
    queries = []; vertices = []; coefficients = []
    for side in (-1, 1):
        up = np.flatnonzero(upper & ((c[:, 2]-midplane)*side > 0))
        low = np.flatnonzero(~upper & ((c[:, 2]-midplane)*side > 0))
        tri = Delaunay(c[low][:, [0, 2]])
        q = c[up][:, [0, 2]]; simplex = tri.find_simplex(q)
        ids = tri.simplices[simplex.clip(0)]
        xy = c[low][ids][:, :, [0, 2]]
        edge = np.max(np.linalg.norm(xy-np.roll(xy, 1, axis=1), axis=2), axis=1)
        valid = (simplex >= 0) & (edge < max_edge)
        transform = tri.transform[simplex[valid]]
        first = np.einsum('nij,nj->ni', transform[:, :2], q[valid]-transform[:, 2])
        weights = np.column_stack([first, 1-first.sum(1)])
        assert weights.min() > -1e-10
        np.testing.assert_allclose(weights.sum(1), 1, atol=1e-12)
        np.testing.assert_allclose(np.sum(c[low][ids[valid]][:, :, [0, 2]]*weights[:, :, None], axis=1), q[valid], atol=1e-12)
        queries.extend(up[valid]); vertices.extend(low[ids[valid]]); coefficients.extend(weights)
    return bi, ii, np.asarray(queries), np.asarray(vertices), np.asarray(coefficients), upper


def jump(pressure, common, queries, vertices, weights, u):
    p = pressure[:, common]
    return 2*((p[:, vertices]*weights[None]).sum(-1)-p[:, queries])/u**2


def integrate_time(value, phase):
    return np.trapezoid(value, phase, axis=0)/(phase[-1]-phase[0])


def analyze(data, u, ref_area, midplane, edge, core, frame, radius, offset, late_start):
    b, i = data['baseline'], data['isolated']
    bi, ii, q, vertices, weights, upper = paired_chart(b, i, midplane, edge)
    phase = i['time']-i['time'][0]
    np.testing.assert_allclose(b['time']-b['time'][0], phase, atol=1e-6)
    values = {name: jump(d['marker_pressure'], idx, q, vertices, weights, u) for (name, d), idx in zip(data.items(), (bi, ii))}
    values['delta'] = values['baseline']-values['isolated']
    xyz = i['centers'][:, ii[q]]; area = np.linalg.norm(i['vectors'][:, ii[q]], axis=-1)
    assert all(np.isfinite(v).all() for v in values.values())
    mean = {k: (v*area).sum(1)/area.sum(1) for k, v in values.items()}
    nearest = int(np.argmin(np.linalg.norm(xyz[frame]-core, axis=1)))
    center = xyz[frame, nearest]
    ref = (np.linalg.norm(xyz[frame]-center, axis=1) <= radius).astype(float)
    target = center+np.asarray(offset)
    control_center = xyz[frame, np.argmin(np.linalg.norm(xyz[frame]-target, axis=1))]
    control = match_area(xyz[frame], area[frame], np.ones(len(q), bool), control_center, (area[frame]*ref).sum())
    late = phase >= late_start
    patches = {}; signals = {}
    for label, mask in [('RV+', ref), ('Upstream', control)]:
        a = area*mask
        signals[label] = {k: np.sum(v*a, axis=1)/a.sum(1) for k, v in values.items()}
        patches[label] = dict(faces=int((mask > 0).sum()), area=float(a[frame].sum()),
            late={k: float(integrate_time(v[late], phase[late])) for k,v in signals[label].items()},
            cycle={k: float(integrate_time(v, phase)) for k,v in signals[label].items()})
    c_b = .7+np.sin(2*np.pi*phase); c_i = -.4+.5*np.cos(2*np.pi*phase)
    shifted = [jump(d['marker_pressure']+.5*u*u*c[:, None], idx, q, vertices, weights, u)
        for d, c, idx in [(b,c_b,bi),(i,c_i,ii)]]
    shift_error = float(np.max(abs((shifted[0]-shifted[1])-values['delta'])))
    np.testing.assert_allclose(shifted[0], values['baseline'], atol=1e-12)
    np.testing.assert_allclose(shifted[1], values['isolated'], atol=1e-12)
    assert shift_error < 1e-12
    gauges = {}
    for name,d in data.items():
        gauges[name] = {}
        for label,mask in [('whole',np.ones(len(d['common_mask']),bool)),('fin',d['pectoral_mask']),('common',d['common_mask'])]:
            net = d['vectors'][:,mask,1].sum(1)/ref_area
            gauges[name][label] = dict(dCL_dConstantCp=net.tolist(), max_abs=float(abs(net).max()))
    total_area = np.linalg.norm(i['vectors'][:,ii[upper]],axis=-1).sum(1)
    result = dict(phase=phase.tolist(), matched_faces=len(ii), paired_upper_faces=len(q),
        upper_area_coverage_range=[float((area.sum(1)/total_area).min()),float((area.sum(1)/total_area).max())],
        cycle_means={k:float(integrate_time(v,phase)) for k,v in mean.items()},
        global_signals={k:v.tolist() for k,v in mean.items()}, patches=patches,
        patch_overlap=float(np.sum(np.minimum(ref,control)*area[frame])/np.sum(ref*area[frame])),
        anchor=center.tolist(),control_center=control_center.tolist(), gauge_shift_max_error=shift_error,gauge_sensitivity=gauges,
        definition='J=(p_lower-p_upper)/(0.5 U^2); DeltaJ=J_body_fin-J_isolated. Lower pressure interpolated in reference x-z material chart; no extrapolation, bounded triangulation edges. Not a direct force integral.')
    return result, values, mean, xyz, area, {'RV+':ref,'Upstream':control}, signals


def plot_result(result, values, mean, xyz, patches, signals, frame, origin, late_start, colors):
    phase=np.array(result['phase']);late=phase>=late_start
    fig,axs=plt.subplots(1,3,figsize=(10.2,3.6))
    fig.subplots_adjust(left=.075,right=.98,bottom=.22,top=.85,wspace=.48)
    for ax,p in zip(axs,'abc'):style_axes(ax,p)
    for name,color in zip(['baseline','isolated'],colors):axs[0].plot(phase,mean[name],color=color,label=name,lw=1.2,zorder=0)
    axs[0].set(xlabel=r'$t/T$',ylabel=r'$\langle C_{p,l}-C_{p,u}\rangle_A$',title='Paired surface pressure difference')
    axs[0].legend(frameon=False,fontsize=8)
    for name,color in zip(patches,colors):axs[1].plot(phase,signals[name]['delta'],color=color,label=name,lw=1.2,zorder=0)
    axs[1].set(xlabel=r'$t/T$',ylabel=r'$\Delta(C_{p,l}-C_{p,u})$',title='Coupling change in pressure jump')
    axs[1].legend(frameon=False,fontsize=8)
    x=xyz[frame]-origin;use=x[:,2]<0
    z=integrate_time(values['delta'][late],phase[late])
    sc=axs[2].scatter(x[use,0],x[use,2],c=z[use],s=6,cmap='RdBu_r',vmin=-.5,vmax=.5,zorder=0)
    for mask in patches.values():
        m=(mask>0)&use;axs[2].scatter(x[m,0],x[m,2],s=10,facecolors='none',edgecolors='k',linewidths=.3,zorder=2)
    axs[2].set(xlabel=r'$x^*$',ylabel=r'$z^*$',title='Late mean jump change',aspect='equal',ylim=(-.235,-.075))
    fig.colorbar(sc,ax=axs[2],shrink=.75,pad=.03,extend='both')
    return fig


def main():
    inputs=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022/loads')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-23-pressure-gauge'
    u,ref_area,midplane=.29,.548382,6.
    core=np.array([10.273124694824219,5.939525127410889,5.912974834442139])
    frame,radius,offset,late_start=12,.015,[-.03,0,0],.75
    edge=.02; sensitivity_edges=[.015,.025]
    origin=np.array([10.,6.,6.]);colors=['#d71920','#087bc1'];dpi=600
    matplotlib.rcParams.update({'font.family':'Times New Roman','font.size':9,'mathtext.fontset':'stix','svg.fonttype':'none','pdf.fonttype':42})
    data,manifest=read_loads(inputs)
    result,*plot_data=analyze(data,u,ref_area,midplane,edge,core,frame,radius,offset,late_start)
    result['edge_sensitivity']={}
    for limit in sensitivity_edges:
        check,*_=analyze(data,u,ref_area,midplane,limit,core,frame,radius,offset,late_start)
        result['edge_sensitivity'][str(limit)]={k:check[k] for k in ['paired_upper_faces','upper_area_coverage_range','cycle_means','patches']}
    result['source_manifest']=manifest
    result['settings']=dict(U=u,reference_area=ref_area,midplane=midplane,max_triangle_edge=edge,frame=frame,radius=radius,offset=offset,late_start=late_start)
    out.mkdir(parents=True,exist_ok=True)
    values,mean,xyz,area,patches,signals=plot_data
    result['figure']=save_hybrid(plot_result(result,values,mean,xyz,patches,signals,frame,origin,late_start,colors),out/'08_pressure_jump',dpi)
    (out/'pressure_jump.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['paired_upper_faces','upper_area_coverage_range','cycle_means','patches','patch_overlap','gauge_shift_max_error']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
