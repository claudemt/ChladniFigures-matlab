# Chladni GUI Studio (restructured)

This version keeps the original plotting titles and figure layout, but reorganizes the codebase more clearly:

- `core/common/` : unified run entry, output discovery, rectangular boundary metadata
- `core/rect/`   : rectangular plate solvers
- `core/circ/`   : circular plate solver
- `app/ui/display/` : shared UI display helpers, color bar, and plot-format helpers
- `app/`         : GUI assembly, controls, preview, export
- `docs/`        : GUI note text

## Rectangular branches

Legacy branches kept as-is in spirit:

- `FFFF` : original sparse free-edge route
- `CCCC` : original clamped finite-difference route
- `SSSS` : Navier exact solution

Analytic Levy-family branches added:

- `SSCC`
- `SSFF`
- `SSSC`
- `SSSF`
- `SSCF`

For these `SS??` cases, the implementation follows the lecture's Levy separation framework with four boundary row operators `RW`, `Rtheta`, `RM`, `RV`, and solves the transcendental determinant condition for each modal family.

## Output folders

Generated figures are written into deterministic parameter folders under `output/`, for example:

- `output/rect_sscc_nu_0.225/`
- `output/circ_free_nu_0.3/`

No timestamp is used. Re-running the same parameter set clears old PNG files in that folder before writing new ones.

## Run

```matlab
main
```
