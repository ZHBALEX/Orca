> v13 修正：按用户确认的 Neumann 压力边界，跨算例单侧压力残差的绝对值不能直接解释为额外吸力；以 [零点无关的上下表面压差分析](../2026-10-23-pressure-gauge/README.md) 为准。旧版开放面片载荷与残差比例不构成独立力归因。

# Draft results and discussion — not submission-ready

Working title: **Spatial organization of body–fin interaction and suction enhancement in an undulating orca configuration**

This is an exploratory English draft based on a single saved cycle and previously selected structures. All evidence has been checked by the assistant against the recorded outputs; accountable author verification and journal-specific review remain pending. No authorship, funding, mesh convergence or global-priority statement is supplied.

## Research question and proposed contribution

How does the presence of an undulating body reorganize pectoral-fin-adjacent vorticity and the spatial distribution of pressure loading? The proposed contribution is a configuration-specific account that combines controlled flow comparisons, a streamwise-vorticity sign-change observation, and spatially resolved pressure accounting. It is not a claim that all body–fin interactions, root vortices or secondary vortices are new.

## Results draft

The existing controlled comparisons identify a positive streamwise branch near the pectoral-fin root and trailing-edge region. Six completed paired releases from identical positions exhibited a negative-to-positive change in the streamwise vorticity component in the coupled configuration, whereas the isolated-fin trajectories retained negative values. A seventh pair did not complete the observation window and was excluded from the endpoint comparison. These observations concern the streamwise component along resolved fluid trajectories, rather than the identity or reversal of an entire material vortex tube. [C01; E01: audit_metrics.json, same_release_pairs]

Across the resolved upper-surface layer, h/L=0.012–0.040, the measured cycle-mean pressure difference between the coupled and isolated-fin configurations was −0.075883. Reconstruction from the velocity-based momentum terms gave −0.074159, corresponding to a 2.27% relative discrepancy in the mean. However, the RMS reconstruction error was 37.94% of the measured fluctuation RMS. The agreement of the cycle means therefore does not establish equally accurate phase-resolved closure. [C02; E01: regions.upper.pressure]

To examine spatial correspondence with the previously selected positive branch, we defined nested material patches around the upper-surface face nearest the reference core at t/T=0.75. Patch membership was fixed throughout the cycle. The smallest patch, r/L=0.015, comprised 28 faces and 2.178% of the valid upper-surface area. During t/T=0.75–1, it accounted for 4.680% of the signed lift increment associated with the incident-pressure-removed diagnostic residual. The corresponding values were 6.090% and 11.934% for r/L=0.025, and 15.089% and 26.893% for r/L=0.040. The result demonstrates concentration of the diagnostic surface-loading residual near the selected branch, but does not isolate a vortex-specific force contribution. [C03; E02: local_pressure_bridge.json, rows]

The spatial correspondence did not imply uniform local momentum balance. In the smallest patch, the late-window incident-removed layer-pressure difference was −0.163445. The Lamb contribution was +0.006158, whereas the kinetic-gradient contribution was −0.195659; the unsteady and viscous contributions were +0.052260 and approximately −0.003088. The reconstructed value was −0.140328, leaving a residual of −0.023118. In the larger r/L=0.040 patch, the Lamb contribution became −0.160611. Thus, the relative sizes and signs of the momentum terms depended on the spatial neighborhood. [C04; E02: rows, late_window_mean]

## Discussion draft

These results motivate a distinction between a flow structure, a spatial pressure association and a force attribution. The coupled configuration supports a distinct streamwise-vorticity evolution near the fin; the neighboring material surface also exhibits concentrated additional suction after removal of the body-only incident-pressure diagnostic. However, the smallest patch does not display the locally Lamb-dominated pressure drop that a simple one-vortex/one-suction-patch picture would predict. The kinetic-gradient and Lamb terms are coupled parts of the momentum equation, not independent causal interventions. Their spatial redistribution must therefore be interpreted together with the unsteady contribution, outer pressure and unresolved near-wall interval.

A stronger mechanism claim requires matched spatial controls and direct comparison of the surface maps of these terms. It also requires adequate native temporal sampling: interpolation substeps cannot replace additional saved flow snapshots. The present single-cycle analysis supports a bounded mechanism candidate and identifies discriminating follow-up tests; it does not establish a new universal vortex class, a body-size effect, or a unique percentage of lift generated by the positive branch.

## Evidence and figure outline

| Claim | Evidence | Intended display | Status |
| --- | --- | --- | --- |
| C01 | E01: same_release_pairs | Current Fig. 2; historical v9/v10 | Checked saved outputs; author review pending |
| C02 | E01: regions.upper | Current Fig. 1 | Newly recomputed |
| C03 | E02: rows | Current Fig. 3 | New exploratory material-patch analysis |
| C04 | E02: rows and central_method_rows | Current Fig. 4 | New local momentum analysis; incomplete closure retained |

Background comparisons to Wang et al. (2019), Liu et al. (2017) and He & Breuer (2026) are documented in the Chinese report with source locations. A final introduction and novelty statement should be drafted only after extending the targeted prior-art comparison and confirming the intended journal. No claim of exhaustive literature coverage is made.
