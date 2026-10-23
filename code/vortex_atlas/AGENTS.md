# Orca postprocessing organization

Follow the user's existing `run__orcaFlukeDiff.py` organization and the
`picar-post-python-style` skill when adding or modifying code here.

- User-facing runnable entry points use `run__*.py` with editable settings at
  the beginning of `main()`: paths, cases, phase/frame, scales, labels, colours,
  limits, and export options. Do not scatter these settings across helpers.
- Separate readers (`read_*.py`), numerical/plot helpers (`fn___*.py`), and
  orchestration. Readers receive paths; plot functions receive arrays and
  options and return figures. Do not import historical analysis runners to
  obtain data or geometry helpers.
- The current three-figure entry point is `run__vortexAtlas.py`. Keep only its
  four required helper modules and these two documentation files here.
- Historical analysis, report generators, compatibility launchers, experiments
  and temporary checks belong outside this repository, under the Orca project's
  `test/body_fin_mechanism_20261022/code_archive/` directory. Do not add them back
  to this folder. The 20261023_before_cleanup archive preserves the earlier code.
- Keep the approved labels TEV2 (legacy data ID W2) and RV+/RV-. Axes show only
  x*, y*, z*, with numeric ticks transformed using the runner's origin and L.
- Scientific plot content, including coloured fields, core/Q contours,
  geometry and colour-bar fills, must be embedded raster images in SVG/PDF.
  Text, ticks, axes, arrows and annotation callouts remain editable vectors.
- A structural refactor must preserve the numerical inputs and rendered
  scientific content. Validate the output against the previous figures and
  inspect SVG layer structure; do not silently change scientific conclusions.
- Keep export packages synchronized with the actual runner and all helpers.
