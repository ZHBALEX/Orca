"""Deterministic sensitivity checks, not tests on independent replicates."""
import numpy as np


def time_mean(value, phase):
    return float(np.trapezoid(value, phase) / (phase[-1] - phase[0]))


def error_metrics(measured, predicted, phase):
    error = predicted - measured
    mean = time_mean(measured, phase)
    rms = np.sqrt(time_mean(error**2, phase))
    fluctuation_rms = np.sqrt(time_mean((measured - mean)**2, phase))
    return dict(measured_mean=mean, predicted_mean=time_mean(predicted, phase),
                signed_mean_error=time_mean(error, phase),
                mean_relative_error=abs(time_mean(error, phase) / mean),
                rms_error=float(rms), fluctuation_rms=float(fluctuation_rms),
                normalized_rms_error=float(rms / fluctuation_rms),
                max_absolute_error=float(np.max(np.abs(error))),
                phase_r_excluding_duplicate_endpoint=float(np.corrcoef(measured[:-1], predicted[:-1])[0, 1]))


def audit(data):
    regions = {}
    for name, row in data['phases'].items():
        phase = row['phase']
        assert len(phase) == 17 and np.allclose(phase, np.linspace(0, 1, 17))
        for kind in ('pressure', 'lift'):
            np.testing.assert_allclose(row[f'{kind}_resolved_prediction'],
                                       sum(row[f'{kind}_{k}'] for k in ('unsteady', 'lamb', 'kinetic', 'viscous')), atol=1e-12)
            np.testing.assert_allclose(row[f'{kind}_surface'],
                                       sum(row[f'{kind}_{k}'] for k in ('outer', 'unresolved_wall', 'resolved_measured')), atol=1e-12)
        quantities = {}
        for key, values in row.items():
            if key == 'phase':
                continue
            mean = time_mean(values, phase)
            quarters = [float(np.trapezoid(values[j:j+5], phase[j:j+5])) for j in (0, 4, 8, 12)]
            np.testing.assert_allclose(sum(quarters), mean, atol=1e-12)
            quantities[key] = dict(trapezoidal_mean=mean,
                                  first_16_mean=float(values[:-1].mean()),
                                  incorrect_17_point_mean=float(values.mean()),
                                  endpoint_difference=float(values[-1]-values[0]),
                                  signed_quarter_integrals=quarters)
        regions[name] = dict(quantities=quantities,
                            pressure=error_metrics(row['pressure_resolved_measured'], row['pressure_resolved_prediction'], phase),
                            lift=error_metrics(row['lift_resolved_measured'], row['lift_resolved_prediction'], phase),
                            lamb_pressure_r=float(np.corrcoef(row['pressure_lamb'][:-1], row['pressure_resolved_measured'][:-1])[0, 1]))
        reference = next(r for r in data['budget']['rows'] if r['region']==name and r['method']=='fourier' and r['h0']==.012 and r['h1']==.04)
        for kind in ('pressure', 'lift'):
            for key, expected in reference[kind]['mean_terms'].items():
                np.testing.assert_allclose(quantities[f'{kind}_{key}']['trapezoidal_mean'], expected, atol=1e-12)
    intervals = []
    for row in data['budget']['rows']:
        if row['region'] != 'upper':
            continue
        p, l = row['pressure']['mean_terms'], row['lift']['mean_terms']
        intervals.append(dict(method=row['method'], h0=row['h0'], h1=row['h1'],
                              area_fraction=row['area_fraction'], mean_pressure=p['resolved_measured'],
                              lamb=p['lamb'], kinetic=p['kinetic'], unsteady=p['unsteady'],
                              mean_relative_error=abs(p['closure_error']/p['resolved_measured']),
                              unfitted_R2=row['pressure']['unfitted_R2'],
                              resolved_lift_fraction=l['resolved_measured']/l['surface'],
                              lamb_lift_fraction=l['lamb']/l['surface']))
    paired = []
    rows = data['rotation']['summary']
    for b in (r for r in rows if r['case']=='baseline'):
        i = next(r for r in rows if r['case']=='isolated' and r['offset']==b['offset'])
        np.testing.assert_allclose(b['initial_position'], i['initial_position'], atol=1e-12)
        paired.append(dict(offset=b['offset'], completed=b['completed'] and i['completed'],
                           baseline_start=b['start_omega'][0], baseline_end=b['end_omega'][0],
                           isolated_start=i['start_omega'][0], isolated_end=i['end_omega'][0]))
    completed = [r for r in paired if r['completed']]
    assert len(completed)==6 and len(paired)==7
    assert all(r['baseline_start']<0<r['baseline_end'] and r['isolated_start']<0 and r['isolated_end']<0 for r in completed)
    return dict(regions=regions, interval_sensitivity=intervals, same_release_pairs=paired,
                scope='Reanalysis of saved derived outputs; no new CFD simulation, independent replicate, or native time/grid convergence test.',
                checks='Phase grid, budget identities, original means, quarter additivity and paired trajectory endpoints passed.')
