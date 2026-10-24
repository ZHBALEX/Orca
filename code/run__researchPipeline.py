"""One command: raw interfaces -> pressure/load chain -> figures -> report.

Local datasets and TecIO dependencies are intentionally not uploaded to GitHub.
Configure local paths in run__rawInterfaceAudit.main and run__integratedChain.main.
"""
from pathlib import Path
import json
import time
import numpy as np
from run__rawInterfaceAudit import main as audit_raw
from run__integratedChain import main as compute_chain
from fn___plotChain import figures
from fn___reportChain import write_report


def main():
    output=Path(__file__).resolve().parents[1]/'reports/2026-10-23-integrated-chain'
    start=time.perf_counter()
    audit_raw()
    compute_chain()
    metrics=json.loads((output/'integrated_chain.json').read_text())
    with np.load(output/'chain_signals.npz') as f:signals={k:f[k] for k in f.files}
    visual=figures(output,metrics,signals)
    write_report(output,metrics)
    (output/'pipeline_run.json').write_text(json.dumps(dict(status='completed',
        elapsed_seconds=time.perf_counter()-start,raw_frames=34,potential_results_audited=102,
        figures=visual,scientific_attribution_closed=False,
        note='Completed analysis pipeline is not equivalent to a closed causal mechanism.'),indent=2),encoding='utf-8')
    print('Completed: raw audit, integrated chain, figures, report.',flush=True)


if __name__=='__main__':main()
