"""Frozen geometry-based controls for the RV+ pressure association."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
from run__localPressureBridge import read_profiles, measure
from fn___spatialEvidence import residual_fields, define_patches, summarize_patches, plot_controls, plot_maps, plot_lamb_split
from fn___plotEvidence import save_hybrid


def main():
    profile_root = Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\reports\orca_pressure_origin_20261022\data')
    output = Path(__file__).resolve().parents[1]/'reports'/'2026-10-23-mechanism-logic'
    core = np.array([10.273124694824219, 5.939525127410889, 5.912974834442139])
    frame, radius, late_start = 12, .015, .75
    offsets = {'x−': [-.03, 0, 0], 'x+': [.03, 0, 0], 'z−': [0, 0, -.03], 'z+': [0, 0, .03]}
    h0, h1, u, period, nu, length = .012, .040, .29, 1., 1/41379.3, 1.
    origin = np.array([10., 6., 6.])
    colors = ['#d71920', '#087bc1', '#a11c68', '#20864c', '#d18917']
    split_colors = ['#087bc1', '#d71920', '#92979b', '#333333']
    map_limit, cmap, dpi = .4, 'RdBu_r', 600
    matplotlib.rcParams.update({'font.family': 'Times New Roman', 'font.size': 10, 'mathtext.fontset': 'stix', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fields, g, manifest = read_profiles(profile_root)
    methods = {m: residual_fields(fields, g, h0, h1, u, period, nu, m) for m in ('fourier', 'central')}
    valid = g['upper'].copy()
    for q in methods.values():
        for v in q.values():
            valid &= np.isfinite(v).all(axis=0)
    patches, metadata = define_patches(g, valid, core, frame, radius, offsets)
    disjoint, disjoint_metadata = define_patches(g, valid, core, frame, radius, offsets, exclude_reference=True)
    result = dict(settings=dict(frame=frame, radius=radius, late_start=late_start, offsets=offsets,
        h0=h0, h1=h1, U=u, period=period, nu=nu, L=length, core_xyz=core.tolist()),
        patches=metadata, source_manifest=manifest, valid_faces=int(valid.sum()), methods={},
        nonoverlap_sensitivity=dict(patches=disjoint_metadata, methods={}))
    for method, q in methods.items():
        rows, signals = summarize_patches(q, g, patches, late_start)
        result['methods'][method] = rows
        result['nonoverlap_sensitivity']['methods'][method] = summarize_patches(q, g, disjoint, late_start)[0]
        old, _, _ = measure(fields, g, core, frame, [radius], h0, h1, u, period, nu, method)
        old_values = old['rows'][1]['late_window_mean']
        for new, previous in [('surface', 'residual'), ('layer', 'layer_residual'), ('prediction', 'predicted_residual'), ('closure', 'closure_residual')]:
            for kind in ('cp', 'cl'):
                np.testing.assert_allclose(rows['RV+']['late'][new+'_'+kind], old_values[previous+'_'+kind], atol=1e-12)
        if method == 'fourier':
            output.mkdir(parents=True, exist_ok=True)
            result['control_figure'] = save_hybrid(plot_controls(g, patches, rows, signals, frame, origin, length, colors), output/'05_spatial_controls', dpi)
            result['map_figure'] = save_hybrid(plot_maps(q, g, valid, frame, late_start, patches, origin, length, map_limit, cmap), output/'06_pressure_maps', dpi)
            result['split_figure'] = save_hybrid(plot_lamb_split(rows, split_colors), output/'07_lamb_midpoint_split', dpi)
    result['checks'] = 'Equal area at anchor phase; surface partition identity; RV+ pressure and force reproduce v11 for both time derivatives; saved Lamb field reproduced from velocity and curl of gradient; midpoint velocity/vorticity split exact on every finite profile.'
    for label, f in patches.items():
        ref_area = np.sum(g['area']*patches['RV+'], axis=1)
        ratio = np.sum(g['area']*f, axis=1)/ref_area
        metadata[label]['area_ratio_range_over_cycle'] = [float(ratio.min()), float(ratio.max())]
    result['limitations'] = 'Exploratory spatial controls; not causal interventions or independent replicates. Patch selection uses geometry only, but RV+ reference and late interval derive from earlier data. Area is exactly matched at phase .75, not enforced at other deformed phases. Closest mesh center may differ from requested offset; overlap reported. No significance test.'
    (output/'spatial_controls.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(dict(patches=metadata, late={k:v['late'] for k,v in result['methods']['fourier'].items()}, checks=result['checks']), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
