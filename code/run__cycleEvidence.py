"""Audit recorded-cycle stability and force sampling; no synthetic CFD fields."""
from pathlib import Path
import hashlib
import json
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import save_hybrid


def read_cycles(path, u, area, steps):
    raw = np.loadtxt(path)
    time = raw[:, 0]
    assert np.all(np.diff(time) > 0)
    idx = np.rint(time * steps).astype(int)
    assert np.max(abs(time - idx / steps)) < 6e-5
    assert np.all(np.diff(idx) == 1)
    cycles = {}
    for start in range(int(np.ceil(time[0])), int(np.floor(time[-1]))):
        mask = (idx >= start * steps) & (idx <= (start + 1) * steps)
        assert mask.sum() == steps + 1
        t = time[mask] - start
        values = raw[mask, 4:7] / (u*u*area)
        cycles[start] = (t, values)
    return cycles, dict(source=path.parent.parent.name + '/forces/' + path.name,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(), rows=len(raw))


def summarize(cycles):
    out = []
    previous = None
    for start, (t, y) in cycles.items():
        row = dict(start=start, mean=np.trapezoid(y, t, axis=0).tolist())
        if previous is not None:
            delta = y - previous
            row['previous_cycle_rms_difference'] = np.sqrt(np.trapezoid(delta**2, t, axis=0)).tolist()
        out.append(row)
        previous = y
    return out


def sampling(t, signal, steps):
    dense = float(np.trapezoid(signal, t))
    out = []
    for intervals in [16, 32, 60, 120, 240, 480]:
        stride = steps // intervals
        ix = np.arange(0, steps + 1, stride)
        # All offsets on a cyclic extension of the half-open recorded cycle.
        # This measures quadrature-phase sensitivity, not extra simulations.
        shifts = [float(signal[offset:steps:stride].mean()) for offset in range(stride)]
        np.testing.assert_allclose(np.mean(shifts), signal[:-1].mean(), atol=1e-14)
        out.append(dict(intervals=intervals, endpoint_trapezoid=float(np.trapezoid(signal[ix], t[ix])),
            cyclic_offset_min=min(shifts), cyclic_offset_max=max(shifts),
            cyclic_offset_mean=float(np.mean(shifts)), dense=dense))
    return out


def main():
    root = Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D')
    out = Path(__file__).resolve().parents[1] / 'reports/2026-10-23-cycle-evidence'
    names = dict(B='Sweep35', I='Sweep35_OnlyFin', N='NoFin', BS='Sweep35_Static', IS='Sweep35_OnlyFin_Static')
    u, area, steps = .29, .548382, 480
    out.mkdir(parents=True, exist_ok=True)
    histories = {}; provenance = {}; stats = {}
    for key, name in names.items():
        histories[key], provenance[key] = read_cycles(root/'Result'/name/'forces/drag_lift_body_001.dat', u, area, steps)
        stats[key] = summarize(histories[key])
    t, b = histories['B'][9]
    i = histories['I'][3][1]; n = histories['N'][9][1]
    signal = b - i - n
    mean = np.trapezoid(signal, t, axis=0)
    np.testing.assert_allclose(np.trapezoid(b, t, axis=0)-np.trapezoid(i, t, axis=0)-np.trapezoid(n, t, axis=0), mean, atol=1e-14)
    late = {key: stats[key][-3:] for key in ['B', 'I', 'N']}
    combinations = np.array([np.array(b['mean'])-i['mean']-np.array(n['mean'])
        for b, i, n in itertools.product(late['B'], late['I'], late['N'])])
    audit = dict(provenance=provenance, cycles=stats,
        columns=['pressure', 'viscous', 'total'], selected_cycle_mean=mean.tolist(),
        late_cycle_cartesian_product=dict(count=len(combinations), minimum=combinations.min(0).tolist(), maximum=combinations.max(0).tolist(),
            interpretation='Sensitivity envelope over last three recorded complete cycles per case; correlated time records, not independent replicates or a confidence interval.'),
        pressure_sampling=sampling(t, signal[:, 0], steps),
        scope='Native whole-system B-I-N force, distinct from common-fin paired pressure; force quadrature cannot certify convergence of field derivatives.')
    (out/'cycle_evidence.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    for key, color in [('B','#0072B2'), ('N','#D55E00')]:
        rows = stats[key]
        axes[0].plot([r['start']+.5 for r in rows], [r['mean'][0] for r in rows], 'o-', label=key, color=color, zorder=0)
    axes[0].set(xlabel='Cycle midpoint, t/T', ylabel='Mean pressure lift coefficient')
    axes[0].legend(frameon=False)
    rows=audit['pressure_sampling']
    axes[1].plot([r['intervals'] for r in rows], [r['endpoint_trapezoid'] for r in rows], 'o-',color='#0072B2',zorder=0)
    axes[1].axhline(mean[0],ls='--',color='black',lw=.8,label='480 intervals',zorder=0)
    axes[1].set(xscale='log',xlabel='Intervals per cycle',ylabel='Mean pressure lift, B - I - N')
    axes[1].set_xticks([16,32,60,120,240,480], ['16','32','60','120','240','480'])
    axes[1].minorticks_off()
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.set_rasterization_zorder(1)
        ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    save_hybrid(fig, out/'16_cycle_sampling', 600)
    plt.close(fig)
    print(json.dumps({k:v for k,v in audit.items() if k not in ['cycles','provenance']},indent=2))


if __name__ == '__main__':
    main()
