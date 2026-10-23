# ORCA-EVIDENCE-012 — spatial specificity and pressure maps

Status: ready within the user's 2026-10-23 instruction to continue Orca research. This records a bounded next analysis; it is not an unbounded autonomous simulation mandate or a scheduled run.

## Objective

Determine whether the extra suction found near the selected RV+ anchor is spatially specific or part of a broader root-region pressure distribution. Explain the change from kinetic-gradient-dominated inner-patch accounting to stronger Lamb contributions over a larger neighborhood.

## Starting evidence

Read README.md, local_pressure_bridge.json and code/run__localPressureBridge.py from the v11 report/commit. Preserve the original selected anchor and the fact that these observations are exploratory. Use the next uploaded commit as the tested baseline; do not invent its SHA before upload.

## Permitted scope

Use existing geometry and saved profile fields; add maintained analysis helpers and figures to Orca. Keep the input data local. No solver settings, remote cluster jobs, manuscript submission or bulk dataset migration.

## Analysis and rival explanations

| Candidate explanation | Discriminating observation | Counter-result / limitation |
| --- | --- | --- |
| Selected branch vicinity has spatially concentrated extra suction | Fixed-size/area neighborhoods centered near RV+ outperform prespecified upstream, downstream and spanwise offsets | Similar or larger concentrations in control neighborhoods weaken specificity |
| General root pressure gradient explains concentration | Concentration varies smoothly with root location even outside branch vicinity | A localized departure remains after geometric matching |
| Neighborhood/mask selection drives the effect | Results change strongly under common all-phase validity or shifted anchor | Stable signs across reasonable masks support robustness, not causality |
| Pressure transport/projection explains local budget change | Same-face maps show kinetic-gradient and Lamb terms occupying different regions, with documented cancellation | Large closure residual prevents a firm local mechanism assignment |

Freeze the offset list and minimum valid area before viewing new patch outcomes; report every attempted neighborhood. Use matched material coordinates, retain geometry and area differences, and avoid selecting the best-looking patch after inspecting pressure. Do not attach p-values to correlated faces or phase samples.

## Acceptance criteria

1. A complete neighborhood table with location, area, missingness, signed residual pressure/force and central/Fourier budgets.
2. One aligned map of surface residual pressure, Lamb, kinetic-gradient and closure residual on identical material faces; one plot comparing the original patch with all controls.
3. Exact surface=outer+layer+wall identities and agreement with v11 for unchanged anchor/mask.
4. A report stating whether spatial specificity strengthened, weakened or remained unresolved, plus an updated English paragraph without inventing a vortex-force fraction.
5. Code, report, small derived tables and report figures uploaded; remote commit verified.

## Deferred checks

Additional native snapshots, independent cycles, grid refinement, body-size or motion-amplitude parameter studies require available simulation outputs or a separate compute setup. They remain open and must not be labeled validated.
