"""Area-matched spatial controls; all patches fixed before reading pressure values."""
import numpy as np
import matplotlib.pyplot as plt
from run__localPressureBridge import normal_terms
from fn___plotEvidence import style_axes


def residual_fields(fields, g, h0, h1, u, period, nu, method):
    ix = np.flatnonzero((g['distances'] >= h0-1e-9) & (g['distances'] <= h1+1e-9))
    b, i, n = [fields[k] for k in ('baseline', 'isolated', 'body_only')]
    q = {'surface': b['surface_cp']-i['surface_cp']-n['surface_cp']}
    p = b['pressure_cp']-i['pressure_cp']-n['pressure_cp']
    q['outer'] = p[:, :, ix[-1]]
    q['layer'] = p[:, :, ix[0]]-q['outer']
    q['wall'] = q['surface']-p[:, :, ix[0]]
    q['lamb'] = -2*np.trapezoid(
        b['lamb_normal'][:, :, ix]-i['lamb_normal'][:, :, ix]-n['lamb_normal'][:, :, ix],
        g['distances'][ix], axis=-1)/u**2
    velocity = [p['velocity'] for p in (b, i, n)]
    omega = []
    for p, v in zip((b, i, n), velocity):
        d = p['grad']
        w = np.stack([d[..., 2, 1]-d[..., 1, 2], d[..., 0, 2]-d[..., 2, 0], d[..., 1, 0]-d[..., 0, 1]], axis=-1)
        omega.append(w)
        rebuilt = np.sum(np.cross(v, w)*g['normal'][:, :, None, :], axis=-1)
        finite = np.isfinite(rebuilt) & np.isfinite(p['lamb_normal'])
        np.testing.assert_allclose(rebuilt[finite], p['lamb_normal'][finite], atol=1e-12)
    ub, ui, un = velocity; wb, wi, wn = omega
    pieces = {'lamb_velocity': np.cross(ub-ui, (wb+wi)/2),
              'lamb_vorticity': np.cross((ub+ui)/2, wb-wi),
              'lamb_body_control': -np.cross(un, wn)}
    for key, value in pieces.items():
        projected = np.sum(value*g['normal'][:, :, None, :], axis=-1)
        q[key] = -2*np.trapezoid(projected[:, :, ix], g['distances'][ix], axis=-1)/u**2
    rebuilt = sum(q[k] for k in pieces)
    finite = np.isfinite(rebuilt) & np.isfinite(q['lamb'])
    np.testing.assert_allclose(rebuilt[finite], q['lamb'][finite], atol=1e-12)
    bt, it, nt = [normal_terms(p, g, ix, u, period, nu, method) for p in (b, i, n)]
    q.update({k: bt[k]-it[k]-nt[k] for k in bt})
    q['prediction'] = sum(q[k] for k in ('lamb', 'kinetic', 'unsteady', 'viscous'))
    q['closure'] = q['layer']-q['prediction']
    return q


def match_area(xyz, area, valid, target_xyz, target_area):
    """Closest valid face first; fractional final-face weight makes area exact."""
    indices = np.flatnonzero(valid)
    order = indices[np.argsort(np.linalg.norm(xyz[indices]-target_xyz, axis=1), kind='stable')]
    cumulative = np.cumsum(area[order])
    stop = int(np.searchsorted(cumulative, target_area))
    if stop >= len(order):
        raise ValueError('Insufficient valid area')
    fraction = np.zeros(len(area))
    fraction[order[:stop]] = 1.
    fraction[order[stop]] = (target_area-(cumulative[stop-1] if stop else 0.))/area[order[stop]]
    np.testing.assert_allclose(np.sum(fraction*area), target_area, rtol=1e-12)
    return fraction


def define_patches(g, valid, core, frame, radius, offsets, exclude_reference=False):
    xyz, area = g['centers'][frame], g['area'][frame]
    anchor = int(np.argmin(np.where(g['upper'], np.linalg.norm(xyz-core, axis=1), np.inf)))
    origin = xyz[anchor]
    reference = valid & (np.linalg.norm(xyz-origin, axis=1) <= radius)
    target_area = area[reference].sum()
    patches = {'RV+': reference.astype(float)}
    centers = {'RV+': origin}
    for label, offset in offsets.items():
        target = origin + np.asarray(offset)
        nearest = int(np.argmin(np.where(valid, np.linalg.norm(xyz-target, axis=1), np.inf)))
        centers[label] = xyz[nearest]
        eligible = valid & ~reference if exclude_reference else valid
        patches[label] = match_area(xyz, area, eligible, centers[label], target_area)
    metadata = {}
    for label, f in patches.items():
        metadata[label] = dict(faces=int(np.count_nonzero(f)), center_xyz=centers[label].tolist(),
            displacement=(centers[label]-origin).tolist(), area=float(np.sum(f*area)),
            reference_overlap_fraction=float(np.sum(np.minimum(f, patches['RV+'])*area)/target_area),
            max_radius=float(np.linalg.norm(xyz[f > 0]-centers[label], axis=1).max()))
    return patches, metadata


def summarize_patches(q, g, patches, late_start):
    phase = g['phase']; late = phase >= late_start
    rows = {}; signals = {}
    for name, fraction in patches.items():
        use = fraction > 0
        a = g['area'][:, use]*fraction[use]
        w = g['weight'][:, use]*fraction[use]
        s = {}
        for key, value in q.items():
            s[key+'_cp'] = np.sum(value[:, use]*a, axis=1)/a.sum(axis=1)
            s[key+'_cl'] = np.sum(value[:, use]*w, axis=1)
        signals[name] = s
        rows[name] = {
            'late': {k: float(np.trapezoid(v[late], phase[late])/(phase[late][-1]-phase[late][0])) for k, v in s.items()},
            'cycle': {k: float(np.trapezoid(v, phase)/(phase[-1]-phase[0])) for k, v in s.items()}}
        np.testing.assert_allclose(s['surface_cp'], s['outer_cp']+s['layer_cp']+s['wall_cp'], atol=1e-12)
    return rows, signals


def plot_controls(g, patches, rows, signals, frame, origin, length, colors):
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.7))
    fig.subplots_adjust(left=.075, right=.99, bottom=.24, top=.85, wspace=.48)
    for ax, p in zip(axs, 'abc'):
        style_axes(ax, p)
    xyz = (g['centers'][frame]-origin)/length
    axs[0].scatter(xyz[g['upper'], 0], xyz[g['upper'], 2], s=2, c='.85', zorder=0)
    for (name, f), color in zip(patches.items(), colors):
        mask = f > 0
        axs[0].scatter(xyz[mask, 0], xyz[mask, 2], s=7, color=color, zorder=0)
        axs[1].plot(g['phase'], signals[name]['surface_cp'], '-', color=color, lw=1.2, label=name, zorder=0)
    axs[0].set(xlabel=r'$x^*$', ylabel=r'$z^*$', title='Matched patches', aspect='equal', ylim=(-.235, -.075))
    axs[1].set(xlabel=r'$t/T$', ylabel=r'$C_{p,B}-C_{p,I}-C_{p,N}$', title='Surface pressure residual')
    axs[1].legend(frameon=False, fontsize=8, ncol=2)
    labels = list(rows)
    axs[2].bar(range(len(rows)), [rows[k]['late']['surface_cp'] for k in labels], color=colors, zorder=0)
    axs[2].set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    axs[2].set(ylabel='Late mean surface residual', title=r'$0.75\leq t/T\leq1$')
    return fig


def plot_maps(q, g, valid, frame, late_start, patches, origin, length, limit, cmap):
    fig, axs = plt.subplots(2, 3, figsize=(9.5, 6))
    fig.subplots_adjust(left=.08, right=.86, bottom=.09, top=.92, wspace=.37, hspace=.40)
    keys = ['surface', 'lamb', 'kinetic', 'unsteady', 'layer', 'closure']
    titles = ['Surface residual', 'Layer: Lamb', 'Layer: kinetic', 'Layer: unsteady', 'Layer: measured', 'Layer: closure error']
    xyz = (g['centers'][frame]-origin)/length
    late = g['phase'] >= late_start; t = g['phase'][late]
    ref = patches['RV+'] > 0
    for ax, panel, key, title in zip(axs.flat, 'abcdef', keys, titles):
        style_axes(ax, panel)
        value = np.trapezoid(q[key][late], t, axis=0)/(t[-1]-t[0])
        sc = ax.scatter(xyz[valid, 0], xyz[valid, 2], c=value[valid], s=8, vmin=-limit, vmax=limit, cmap=cmap, zorder=0)
        ax.scatter(xyz[ref, 0], xyz[ref, 2], s=13, facecolors='none', edgecolors='k', linewidths=.25, zorder=2)
        ax.set(xlabel=r'$x^*$', ylabel=r'$z^*$', title=title, aspect='equal', ylim=(-.235, -.075))
    cb = fig.colorbar(sc, cax=fig.add_axes([.90, .20, .018, .58]), extend='both')
    cb.set_label('Late mean pressure coefficient contribution')
    return fig


def plot_lamb_split(rows, colors):
    fig, ax = plt.subplots(figsize=(7.3, 3.7))
    fig.subplots_adjust(left=.12, right=.98, bottom=.20, top=.80)
    style_axes(ax, 'a')
    ax.texts[0].set_position((-.08, 1.06))
    x = np.arange(len(rows))
    for shift, key, label, color in zip(
            [-.27, -.09, .09, .27],
            ['lamb_velocity_cp', 'lamb_vorticity_cp', 'lamb_body_control_cp', 'lamb_cp'],
            ['Velocity difference', 'Vorticity difference', 'Body control', 'Sum'], colors):
        ax.bar(x+shift, [r['late'][key] for r in rows.values()], width=.18, label=label, color=color, zorder=0)
    ax.set_xticks(x, list(rows))
    ax.set(ylabel='Late mean layer-pressure contribution', title='Exact midpoint expansion of the Lamb residual')
    ax.legend(frameon=False, ncol=2, fontsize=8, loc='lower right')
    return fig
