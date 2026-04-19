# ChladniFigures-matlab

A MATLAB project for generating **Chladni figures** (nodal-line patterns of vibrating plates). The code is built around the **Kirchhoff--Love thin-plate model**, computes plate eigenmodes for several geometries and boundary conditions, and exports the zero level sets of those modes as PNG images.

The repository is useful in three ways:

- as a GUI tool for quickly generating rectangular, circular, and annular Chladni patterns;
- as a compact teaching/research code for thin-plate eigenvalue problems;
- as a starting point for extending solvers, geometries, or boundary-condition families.

---

## 1. What the project does

### Supported domains

- **Rectangle**: `rect`
- **Solid disk**: `circ`
- **Annulus**: `annulus`

### Supported boundary conditions

#### Rectangular plate
Rectangular boundaries are encoded in **ULDR** order:

- `U` = up (top)
- `L` = left
- `D` = down (bottom)
- `R` = right

Each side can be:

- `C` = clamped
- `S` = simply supported
- `F` = free

Typical rectangular cases handled in the code include:

- `FFFF`
- `CCCC`
- `SSSS`
- `SSCC`
- `SSFF`
- `SSSC`
- `SSSF`
- `SSCF`

In practice:

- `SSSS` uses a Navier sine-series solution;
- `SS??` families are treated by a Levy-type semi-analytical method;
- more general rectangular cases are handled by Ritz / discrete eigen-solvers.

#### Disk / annulus
- Solid disks support `C`, `S`, and `F`.
- Annuli use a two-letter **outer-inner** code, for example:
  - `CC`, `CS`, `CF`
  - `SC`, `SS`, `SF`
  - `FC`, `FS`, `FF`

---

## 2. How to run it

Start the GUI with:

```matlab
main
```

or

```matlab
launch_chladni_studio(pwd)
```

### Main GUI parameters

- `domain`: `rect / circ / annulus`
- `boundary`: boundary-condition code
- `nu`: Poisson ratio \(\nu\)
- `number of modes`: number of modes to export
- `grid size`: plot resolution
- `xi_0`: geometry ratio parameter

The meaning of `xi_0` depends on the domain:

- for **rect**: with fixed \(a=2\), the code sets \(b = 2\xi_0\), so \(\xi_0 = b/a\);
- for **annulus**: \(\xi_0 = R_0/R\);
- for **circ**: this parameter is disabled and internally taken as \(0\).

### Output

- temporary previews are written to `.cache/`
- exported images are copied to `output/`
- filenames encode geometry, boundary condition, Poisson ratio, and `xi_0`

Examples:

```text
rect-FFFF-nu0.225-xi0.45-mode1.png
annulus-CF-nu0.3-xi0.4-mode2,1.png
```

---

## 3. Governing equations

The project uses the **Kirchhoff--Love thin-plate equation** for a homogeneous, isotropic plate:

\[
\rho h\, \frac{\partial^2 w}{\partial t^2} + D \nabla^4 w = q(x,y,t),
\qquad
D = \frac{E h^3}{12(1-\nu^2)}.
\]

Here:

- \(w(x,y,t)\) is the transverse displacement,
- \(\rho\) is density,
- \(h\) is thickness,
- \(E\) is Young's modulus,
- \(\nu\) is Poisson's ratio,
- \(D\) is the bending stiffness,
- \(\nabla^4 = (\nabla^2)^2\) is the biharmonic operator.

For free vibration, write

\[
w(x,y,t)=W(x,y)e^{i\omega t},
\]

which gives the eigenvalue problem

\[
D\nabla^4 W = \rho h\,\omega^2 W.
\]

The code therefore computes:

- the eigenfunctions \(W\) (mode shapes), and
- the eigenvalues proportional to \(\omega^2\).

A Chladni figure is obtained from the **nodal set**

\[
W(x,y)=0.
\]

Those zero contours are what the exported images display.

### Boundary conditions in continuous form

Let \(n\) be the outward normal direction and \(\tau\) the tangential direction. Then the classical plate boundary conditions are:

#### Clamped (`C`)
\[
W = 0,
\qquad
\frac{\partial W}{\partial n}=0.
\]

#### Simply supported (`S`)
\[
W = 0,
\qquad
M_{nn}=0.
\]

#### Free (`F`)
\[
M_{nn}=0,
\qquad
V_n=0.
\]

Here \(M_{nn}\) is the normal bending moment and \(V_n\) is the effective transverse shear force.

---

## 4. Algorithms used in the project

This repository does **not** use a single solver for every case. Instead, it switches to a method that fits the geometry and boundary-condition family.

### 4.1 Rectangle: Navier solution for `SSSS`

For a plate simply supported on all four sides, the code calls:

```matlab
solve_rect_navier_ssss
```

The modal basis is the classical double-sine family:

\[
W_{mn}(x,y)=\sin\!\left(\frac{m\pi x}{a}\right)
\sin\!\left(\frac{n\pi y}{b}\right),
\]

with spectral parameter ordered from

\[
\lambda_{mn}=\left(\frac{m\pi}{a}\right)^2+\left(\frac{n\pi}{b}\right)^2.
\]

This branch is the cleanest and fastest because the separated modes are explicit.

### 4.2 Rectangle: Levy-type method for `SS??` families

For cases with two opposite simply supported edges and mixed conditions on the other pair, the code uses:

```matlab
solve_rect_levy_family
```

The idea is:

1. expand in sines along the simply supported direction;
2. represent the other direction by combinations of trigonometric and hyperbolic functions;
3. impose the remaining edge conditions through a small characteristic matrix;
4. locate eigenvalues from zeros of the resulting determinant.

In the implementation, roots are detected by scanning the characteristic quantity, identifying sign changes or local minima, and then refining them numerically using `fzero` / `fminbnd`.

### 4.3 Rectangle: Ritz / discrete solvers for general cases

For more general rectangular boundaries, the repository includes:

```matlab
solve_rect_ritz_general
solve_rect_free_ritz_general
solve_rect_free_sparse
solve_rect_clamped_fd_highres
```

These branches assemble approximate stiffness/mass operators and solve matrix eigenvalue problems. In broad terms:

- **Ritz** branches use trial functions chosen to satisfy or approximate the required edge behavior;
- **sparse / finite-difference** branches discretize the plate equation directly on a grid and solve a large sparse eigenproblem.

These solvers are less explicit than the Navier or Levy cases, but they make mixed and fully free/clamped rectangular problems tractable.

### 4.4 Disk and annulus: analytic radial basis + numerical root search

Circular and annular cases are handled in:

```matlab
chladni_formula_circ
```

The method uses separation in polar coordinates. For angular order \(m\), the radial dependence is built from Bessel / modified-Bessel families (depending on the factorization of the biharmonic operator), and the boundary conditions lead to a small analytic system.

For annuli, the code explicitly documents a **two-stage analytic root-finding strategy**:

1. coarse candidate detection using determinant sign changes and minima of the smallest singular value of the boundary matrix;
2. high-accuracy polishing of each root with `fzero` or `fminbnd`.

This makes the circular solver much more stable than naive determinant sampling alone.

---

## 5. Physical meaning of the pictures

Each exported image is a scalar field plot of a mode shape together with its nodal contour \(W=0\). The dark contour lines are the nodal lines where the plate does not move in that mode.

In a physical Chladni experiment, sand accumulates near those nodal lines because the vibration amplitude is smallest there.

---

## 6. Repository structure

```text
main.m                         GUI entry point
app/                           GUI construction and callbacks
core/rect/                     rectangular-plate solvers
core/circ/                     disk / annulus solver
core/common/                   shared utilities and boundary metadata
docs/plate_vibration.tex       extended notes on plate vibration theory
.cache/                        temporary preview images
output/                        exported images
```

---

## 7. Notes and limitations

- The project focuses on **mode generation and visualization**, not on full experimental identification or forced-response simulation.
- Different boundary-condition families are solved by different methods, so numerical behavior and mode ordering may vary slightly across domains.
- For circular problems, the displayed annular and disk modes are real-valued representatives of separated eigenmodes.
- The code is best viewed as a compact computational project connecting classical thin-plate theory to visible nodal patterns.

---

## 8. Suggested use

This repository is a good fit for:

- classroom demonstrations of thin-plate eigenmodes,
- quick generation of Chladni-style figures for reports or talks,
- testing how geometry, aspect ratio, Poisson ratio, and boundary conditions affect nodal lines,
- extending the code toward other geometries or numerical plate solvers.
