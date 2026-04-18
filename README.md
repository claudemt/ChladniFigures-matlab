# ChladniFigures-matlab

MATLAB project for generating Chladni figures of vibrating plates.

## Main formulas

The code solves the thin-plate vibration eigenproblem
\[
D\nabla^4 w = \rho h\,\omega^2 w,
\]
and plots nodal patterns from the mode shape \(w\).

For rectangular plates, the dimensionless eigenvalue shown in figure titles is
\[
\Lambda \propto \omega^2,
\]
with boundary-dependent modal equations.

For circular plates, the radial part uses Bessel / modified Bessel functions,
typically of the form
\[
J_m(\beta r) + C I_m(\beta r),
\]
and eigenvalues are found from the boundary characteristic equation for:

- clamped
- simply supported
- free

## Boundary conventions

Rectangular branches:

- `FFFF`: free on all four edges
- `CCCC`: clamped on all four edges
- `SSSS`: simply supported on all four edges
- `SSCC`, `SSFF`, `SSSC`, `SSSF`, `SSCF`: Levy-type mixed-edge cases

Circular branches:

- `free`
- `simply`
- `clamped`

## Basic conventions

- `nu` = Poisson ratio
- `k` = number of modes to export
- `n` = grid resolution for plotting
- figures are written as PNG files
- nodal lines are the zero contour of the computed mode shape
- repeated runs overwrite/refresh output in the corresponding output folder

## Run

Launch the app:
```matlab
main
```

Or call the core solvers directly after adding the project to the MATLAB path.

## Folder notes

- `app/`: GUI
- `core/rect/`: rectangular plate solvers
- `core/circ/`: circular plate solver
- `core/common/`: shared helpers
- `docs/`: notes and references

## Tips

- Increase `n` for smoother nodal lines.
- Increase `k` to export more modes.
- If a circular case misses modes, enlarge the search range in the circular solver.
