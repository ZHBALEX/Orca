# Orca GitHub workflow

- Project: Pj_Orca_3D.
- Repository: https://github.com/ZHBALEX/Orca (user-created).
- Base branch: main.
- User instruction, 2026-10-23: keep Orca in its own repository and synchronize agent-written code and updates to GitHub after each completed task.
- Upload scope (user clarification, 2026-10-23): only code authored or modified by the assistant for Orca, and corresponding reports for remote review. Keep minimal repository/workflow documentation. Do not bulk migrate existing user code, research documents, input datasets, simulation outputs, or standalone figures. Include figures only as necessary evidence within a report.
- Initial onboarding commit may use main; subsequent work uses codex/ task branches and pull requests. Do not merge without authorization.
- Task authority: direct user requests. No issue or NEXT.md has yet been designated as an executable task source.
- No scheduled/background synchronization was requested.
- Build/test commands: not yet established; verify commands for each scoped code change. Repository setup is not a CFD validation.
- Repository created and cloned by the user on 2026-10-23. Source code checkout: GitHub/Orca; research data remains in the existing Pj_Orca_3D directory. Existing research scripts have not been bulk migrated.

- Research continuation authorized on 2026-10-23: organize existing assistant reports, analyze the evidence and novelty, create scientific figures and an English draft, and progress the Orca mechanism study.
- Current research entry: reports/LATEST.md. The bounded next analysis is documented in reports/2026-10-23-evidence-audit/NEXT.md.
- Validation commands: python code/run__evidenceAudit.py; python code/run__localPressureBridge.py. Both passed using the local scientific environment; input profile data remains local.
- code/vortex_atlas is a copied maintained plotting snapshot; its pipeline was not rerun during the v11 audit.

- v12 completed ORCA-EVIDENCE-012: literature-informed spatial controls and exact midpoint Lamb expansion. Latest report: reports/2026-10-23-mechanism-logic/README.md. New validation: python code/run__spatialControls.py; local consistency checks passed. ORCA-EVIDENCE-013 remains a defined next analysis, not a completed result.

- User clarification, 2026-10-23: pressure has Neumann boundary conditions. Prioritize within-case upper/lower pressure differences or explicitly referenced relative pressure, then compare cases. Do not interpret cross-case single-surface absolute pressure residuals or open-surface pressure integrals as gauge-independent load changes. v13 implements material-chart pressure jumps and independent time-dependent gauge-shift tests; full junction/load closure remains pending.

- v14: run__researchPipeline.py completed raw B/I interfaces (34 phases), dense closed-system loads, gauge-consistent regional pressure/load accounting, paired normal-momentum terms, geometry checks and audits of 102 existing potentials plus 98 transport windows. Report: reports/2026-10-23-integrated-chain/README.md. Remaining physical attribution is explicitly unclosed; no new CFD simulation or elliptic solve claimed.

- v15 completed: run__liftCauseControls.py compares four existing moving/static body/isolated cases using Neumann-invariant material pressure pairs, checks mean shape, pressure/projection covariance, span/phase contributions, 240 section pressure profiles and 102 saved source-term records. Latest explanation is reports/2026-10-23-lift-cause-controls/README.md. No new CFD runs; single-vortex attribution remains distinct from motion-dependent regional pressure enhancement.
