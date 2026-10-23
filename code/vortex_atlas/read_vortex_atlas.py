"""Read existing vortex anchors, field caches and surface meshes; no analysis rerun."""
import json
from pathlib import Path
import numpy as np


def read_vortex_atlas(identity_file, volume_files, geometry_files, phase):
    rows = json.loads(Path(identity_file).read_text(encoding='utf-8'))['rows']
    if any(not np.isclose(row['phase'], phase) for row in rows):
        raise ValueError('Anchor phase differs from requested phase; select matching identities and fields.')
    volumes, geometry = {}, {}
    for case, path in volume_files.items():
        with np.load(path, allow_pickle=False) as data:
            volumes[case] = {key: data[key] for key in ['x', 'y', 'z', 'omega', 'Q', 'valid', 'fluid']}
    for case, path in geometry_files.items():
        with np.load(path, allow_pickle=False) as data:
            geometry[case] = {key: data[key] for key in ['faces', 'pectoral_mask', 'points_history']}
    return rows, volumes, geometry
