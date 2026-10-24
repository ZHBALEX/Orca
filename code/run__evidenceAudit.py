"""Reproduce the v11 evidence audit from local saved Orca outputs."""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
from read_evidence import read_study
from fn___auditEvidence import audit
from fn___plotEvidence import plot_budget, plot_branch, save_hybrid


def main():
    study_root = Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    output = Path(__file__).resolve().parents[1] / 'reports' / '2026-10-23-evidence-audit'
    dpi = 600
    velocity_scale, length_scale = .29, 1.0  # Original study config: U, LENGTH.
    palette = dict(red='#d71920',blue='#087bc1',grey='#70777f')
    matplotlib.rcParams.update({'font.family':'Times New Roman','font.size':10,'mathtext.fontset':'stix',
                               'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7})
    output.mkdir(parents=True,exist_ok=True)
    data,manifest=read_study(study_root)
    results=audit(data)
    results['runtime']=dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__)
    results['source_manifest']=manifest
    results['figure_structure']={}
    for name,figure in [('01_pressure_evidence',plot_budget(data,results,palette)),('02_branch_evidence',plot_branch(data,results,palette,velocity_scale,length_scale))]:
        results['figure_structure'][name]=save_hybrid(figure,output/name,dpi)
    (output/'audit_metrics.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:results['regions']['upper'][k] for k in ['pressure','lift','lamb_pressure_r']},indent=2))
    print('All audit checks passed. Figures and metrics:',output)


if __name__=='__main__':
    main()
