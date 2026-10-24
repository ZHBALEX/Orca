"""Re-read every saved B/I pressure interface; keep bulky audit cache local."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np


def audit_case(folder, data, api, grid_reader, force_reader, pv):
    grids=grid_reader(folder,api)
    cells=np.column_stack([np.full(len(data['faces']),3),data['faces']]).ravel()
    rows=[]
    for j,step in enumerate(data['steps']):
        flux=force_reader(folder,int(step),api,grids)
        mesh=pv.PolyData(data['points_history'][j],cells)
        idx=np.asarray(mesh.find_closest_cell(flux['centers']))
        mapped=np.zeros_like(data['native_force'][j]);np.add.at(mapped,idx,2*flux['force'])
        error=float(np.max(abs(mapped-data['native_force'][j])))
        assert error<1e-10, ('raw/cache mismatch',step,error)
        direction=np.sign(flux['solid_index']-flux['fluid_index'])
        signed_y=flux['area']*direction[:,1]
        row=dict(step=int(step),cache_error=error,regions={})
        for name,mask in [('whole',np.ones(len(data['faces']),bool)),('fin',data['pectoral_mask']),('common',data['common_mask'])]:
            select=mask[idx]
            row['regions'][name]=dict(pressure_force_y=float(2*flux['force'][select,1].sum()),
                signed_area_y=float(signed_y[select].sum()),interfaces=int(select.sum()))
        rows.append(row)
        print(folder.name,int(step),'cache match',f'{error:.2e}',flush=True)
    return rows


def main():
    post=Path(r'C:\Users\ZHB\Documents\GitHub\Picar_Post')
    project=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D')
    loads=project/'test/body_fin_mechanism_20261022/loads'
    output=Path(__file__).resolve().parents[1]/'reports/2026-10-23-integrated-chain'
    cases={'baseline':'Sweep35','isolated':'Sweep35_OnlyFin'}
    sys.path.insert(0,str(post))
    import pyvista as pv
    from tecio_read import TecIO
    from marker_separation.orca.read_cartesian import initialize_grids,read_pressure_faces
    api=TecIO();result={}
    output.mkdir(parents=True,exist_ok=True)
    for name,case in cases.items():
        with np.load(loads/f'{name}.npz') as f:data={k:f[k] for k in f.files}
        result[name]=audit_case(project/'Result'/case,data,api,initialize_grids,read_pressure_faces,pv)
        (output/'raw_interface_audit.partial.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    result['reader_provenance']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [post/'tecio_read.py',post/'marker_separation/orca/read_cartesian.py']}
    result['limitation']='Reconstructs exported P/iblank/AMR_status interfaces. Six internal boundary masks, body ownership and solver double-precision metrics were not exported; agreement with mapped cache is not bitwise identity to solver internals.'
    output.mkdir(parents=True,exist_ok=True)
    (output/'raw_interface_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (output/'raw_interface_audit.partial.json').unlink(missing_ok=True)


if __name__=='__main__':main()
