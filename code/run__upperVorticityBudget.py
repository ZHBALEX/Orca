"""Eulerian vorticity balance at matched moving upper-front sampling points.

Temporal derivatives use neighbouring native frames at FIXED physical points.
Projection onto the instantaneous local span is applied after vector operations.
This is not the derivative of a moving-window circulation.
"""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import binary_erosion


def main():
    study=Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    sys.path.insert(0,str(study/'code_archive/20261023_before_cleanup'))
    from read_volume import fluid_gradient
    rows=[];nu=1/41379.3
    for case in ['baseline','isolated']:
        for phase in range(1,16):
            with np.load(study/'vortex_source_v8/volumes'/f'{case}_{phase:02d}.npz') as f:d=dict(f)
            axes=(d['z'],d['y'],d['x']);fluid=d['fluid'].astype(bool)
            gu,vu=fluid_gradient(d['velocity'],axes,fluid)
            go,vo=fluid_gradient(d['omega'],axes,d['valid'].astype(bool))
            sec=dict(np.load(study/'vortex_structures_v4'/f'baseline_{phase:02d}_sections.npz'))
            queries=[];directions=[]
            for j in range(3):
                key=f'{j}_1_';xx,hh=np.meshgrid(np.linspace(.05,.5,15),np.linspace(.016,.04,8))
                x=sec[key+'leading_x']+xx*sec[key+'chord'];y=np.interp(x,sec[key+'x'],sec[key+'upper'])+hh
                local=np.stack([x,y,np.zeros_like(x)],axis=-1).reshape(-1,3)
                queries.append(sec[key+'origin']+local@sec[key+'basis'].T);directions.append(sec[key+'basis'][:,2])
            points=np.concatenate(queries)
            def sample(a):return RegularGridInterpolator(axes,a,bounds_error=False,fill_value=np.nan)(points[:,::-1])
            valid=sample((vu&vo&binary_erosion(d['valid'].astype(bool),iterations=2)).astype(float))>.999
            v=sample(d['velocity']);om=sample(d['omega']);g=sample(gu);og=sample(go)
            advection=-np.einsum('nj,nij->nji',v,og)
            stretching=np.einsum('nj,nij->nji',om,g)
            divergence=-om*np.trace(g,axis1=1,axis2=2)[:,None]
            lap=np.zeros_like(d['omega'])
            for a in range(3):lap+=np.gradient(np.gradient(d['omega'],axes[a],axis=a,edge_order=2),axes[a],axis=a,edge_order=2)
            diffusion=nu*sample(lap)
            neighbours=[]
            for frame in [phase-1,phase+1]:
                with np.load(study/'vortex_source_v8/volumes'/f'{case}_{frame:02d}.npz') as n:
                    grid=(n['z'],n['y'],n['x'])
                    neighbours.append(RegularGridInterpolator(grid,n['omega'],bounds_error=False,fill_value=np.nan)(points[:,::-1]))
                    valid &= RegularGridInterpolator(grid,n['valid'].astype(float),bounds_error=False,fill_value=0)(points[:,::-1])>.999
            derivative=(neighbours[1]-neighbours[0])/(2/16)
            for j,direction in enumerate(directions):
                sl=slice(j*120,(j+1)*120);ok=valid[sl]
                projected=dict(time=derivative[sl]@direction,divergence=divergence[sl]@direction,conditional_diffusion=diffusion[sl]@direction)
                for k,axis in enumerate('xyz'):
                    projected['advection_'+axis]=advection[sl,k]@direction
                    projected['stretch_tilt_'+axis]=stretching[sl,k]@direction
                projected['residual']=projected['time']-sum(a for key,a in projected.items() if key!='time')
                rows.append(dict(case=case,phase=phase/16,span=[.2,.3,.4][j],valid_fraction=float(ok.mean()),
                    means={key:float(a[ok].mean()) if ok.any() else None for key,a in projected.items()},
                    rms={key:float(np.sqrt(np.mean(a[ok]**2))) if ok.any() else None for key,a in projected.items()}))
            print(case,phase,'done',flush=True)
            del gu,go,lap,d
    (out/'upper_vorticity_budget.json').write_text(json.dumps(dict(rows=rows,equation='partial_t omega = -u_j partial_j omega + omega_j partial_j u - omega div(u) + nu laplacian(omega)',
        scope='Matched baseline-section coordinates, x/c=.05–.5, h=.016–.04 L. Interior native phases only; no assumed periodic endpoints. Each row uses valid points in that case, so cross-case differences require matching coverage. Gradients sampled then multiplied; residual includes interpolation and discretization errors. Cartesian contributions projected on instantaneous span. These are local Eulerian tendencies, not integrated circulation sources or force fractions.',
        viscosity=nu,viscosity_provenance='Conditional on available input.dat Re=41379.3; run-copy identity not independently archived.'),indent=2,allow_nan=False),encoding='utf-8')


if __name__=='__main__':main()
