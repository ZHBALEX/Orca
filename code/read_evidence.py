"""Read the saved Orca evidence; never change the original analyses."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def read_phase_table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    grouped = {}
    for region in sorted({row['region'] for row in rows}):
        subset = sorted((r for r in rows if r['region'] == region), key=lambda r: float(r['phase']))
        grouped[region] = {key: np.array([float(r[key]) for r in subset])
                           for key in subset[0] if key != 'region'}
    return grouped


def read_study(root):
    root = Path(root)
    files = {
        'budget': 'pressure_origin_v7/data/paired_pressure_budget.json',
        'phases': 'pressure_origin_v7/data/paired_pressure_budget_phases.csv',
        'incident': 'pressure_origin_v7/data/incident_pressure_phases.csv',
        'rotation': 'vortex_discovery_v9/rotation_reversal_verification.json',
        'core': 'vortex_discovery_v9/report_metrics.json',
    }
    data = {key: (read_phase_table(root / rel) if rel.endswith('.csv') else read_json(root / rel))
            for key, rel in files.items()}
    manifest = [{'id': key, 'path_relative_to_study': rel,
                 'sha256': hashlib.sha256((root / rel).read_bytes()).hexdigest()}
                for key, rel in files.items()]
    return data, manifest
