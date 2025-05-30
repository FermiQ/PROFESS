# Overview

The `NMS` module is a Fortran implementation of the PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) package, primarily developed by Fred Fritsch. This package provides robust routines for constructing and evaluating one-dimensional interpolants that can preserve monotonicity if present in the input data. It can also be used for standard cubic spline interpolation with various boundary conditions.

The core idea is to represent the interpolating function as a set of cubic polynomials, each defined on an interval between data points. The PCHIP method focuses on setting the derivative values at each data point in such a way that the resulting piecewise cubic function is monotone if the data itself is monotone. For spline interpolation, it determines the derivatives to ensure higher continuity (continuous first and second derivatives).

**Key public routines are `pchez` (to set up the interpolant by computing derivatives) and `pchev` (to evaluate the interpolant and its derivative).**

**References:**
- Fred Fritsch, Ralph Carlson, "Monotone Piecewise Cubic Interpolation," SIAM Journal on Numerical Analysis, Vol. 17, No. 2, April 1980, pp. 238-246.
- Fred Fritsch, Judy Butland, "A Method for Constructing Local Monotone Piecewise Cubic Interpolants," SIAM Journal on Scientific and Statistical Computing, Vol. 5, No. 2, 1984, pp. 300-304.
- Carl de Boor, "A Practical Guide to Splines," Springer-Verlag, 1978.

# Key Components

- **`MODULE NMS`**: The main container for the PCHIP routines.

- **`SUBROUTINE pchez(n, x, f, d, spline, wk, lwk, ierr)`**:
  - **Description:** An "easy-to-use" driver for setting up either a piecewise cubic Hermite interpolant (if `spline` is `.FALSE.`) or a cubic spline interpolant (if `spline` is `.TRUE.`). It computes the derivative values `d(i)` at each data point `x(i)`.
  - **Arguments:**
    - `n :: INTEGER, INTENT(IN)`: Number of data points (must be >= 2).
    - `x(n) :: REAL(KIND=DP), INTENT(IN)`: Strictly increasing independent variable values.
    - `f(n) :: REAL(KIND=DP), INTENT(IN)`: Function values `f(i) = F(x(i))`.
    - `d(n) :: REAL(KIND=DP), INTENT(OUT)`: Calculated derivative values at `x(i)`.
    - `spline :: LOGICAL, INTENT(IN)`: If `.TRUE.`, computes derivatives for a cubic spline (using `pchsp`). If `.FALSE.`, computes derivatives for a PCHIP (monotone if data is monotone, using `pchim`).
    - `wk(lwk) :: REAL(KIND=DP), INTENT(OUT)`: Workspace array, needed only if `spline` is `.TRUE.`.
    - `lwk :: INTEGER, INTENT(IN)`: Length of `wk` (must be at least `2*n` if `spline` is `.TRUE.`).
    - `ierr :: INTEGER, INTENT(OUT)`: Error flag.

- **`SUBROUTINE pchev(n, x, f, d, nval, xval, fval, dval, ierr)`**:
  - **Description:** Evaluates the piecewise cubic function (defined by `x`, `f`, and pre-computed derivatives `d` from `pchez` or other PCHIP setup routines) and its first derivative at an array of specified points `xval`.
  - **Arguments:**
    - `n :: INTEGER, INTENT(IN)`: Number of data points.
    - `x(n), f(n), d(n) :: REAL(KIND=DP), INTENT(IN)`: Original data points and their derivatives.
    - `nval :: INTEGER, INTENT(IN)`: Number of points at which to evaluate.
    - `xval(nval) :: REAL(KIND=DP), INTENT(IN)`: Points at which to evaluate the interpolant.
    - `fval(nval) :: REAL(KIND=DP), INTENT(OUT)`: Interpolated function values at `xval`.
    - `dval(nval) :: REAL(KIND=DP), INTENT(OUT)`: Interpolated derivative values at `xval`.
    - `ierr :: INTEGER, INTENT(OUT)`: Error flag (positive value indicates number of extrapolation points).

- **`SUBROUTINE pchim(n, x, f, d, incfd, ierr)`**:
  - **Description:** Sets derivative values `d` for a piecewise cubic Hermite interpolant, specifically aiming to preserve monotonicity if present in the data `(x,f)`. It handles boundary conditions using local quadratic formulas unless monotonicity is affected.
  - **Arguments:** `n`, `x`, `f`, `d`, `incfd` (increment for f and d, for multi-dimensional array layout), `ierr`.

- **`SUBROUTINE pchsp(ic, vc, n, x, f, d, incfd, wk, nwk, ierr)`**:
  - **Description:** Sets derivative values `d` for a cubic spline interpolant with specified boundary conditions. It solves a tridiagonal linear system for the derivative values.
  - **Arguments:**
    - `ic(2) :: INTEGER, INTENT(IN)`: Specifies boundary conditions at `x(1)` and `x(n)`. (0: "not-a-knot", 1: given derivative, 2: given second derivative, 3: 3-point diff, 4: 4-point diff).
    - `vc(2) :: REAL(KIND=DP), INTENT(IN)`: Values for boundary conditions if `ic(1)` or `ic(2)` is 1 or 2.
    - `n`, `x`, `f`, `d`, `incfd`, `wk` (workspace, size `nwk`), `nwk`, `ierr`.

- **`SUBROUTINE pchfd(n, x, f, d, incfd, skip, ne, xe, fe, de, ierr)`**:
  - **Description:** Evaluates a piecewise cubic Hermite function (defined by `x,f,d`) and its first derivative at an array of points `xe`. Called by `pchev`. Handles extrapolation if `xe` points are outside `[x(1), x(n)]`.
  - **Arguments:** `n`, `x`, `f`, `d`, `incfd`, `skip` (logical for error checks), `ne`, `xe`, `fe` (output func values), `de` (output deriv values), `ierr`.

- **`SUBROUTINE chfdv(x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr)`**:
  - **Description:** Evaluates a single cubic polynomial (defined by Hermite data `f1,d1` at `x1` and `f2,d2` at `x2`) and its derivative at an array of points `xe`.
  - **Arguments:** Endpoint data `x1,x2,f1,f2,d1,d2`, number of evaluation points `ne`, evaluation points `xe`, output `fe, de`, `next(2)` (count of extrapolation points), `ierr`.

- **`FUNCTION pchdf(k, x, s, ierr) RESULT(derivative_approx)`** (Internal to `pchsp`):
  - **Description:** Approximates a derivative at `x(k)` using `k`-point divided differences of slopes `s`. Used by `pchsp` for certain boundary condition calculations.
  - **Arguments:** `k` (order), `x(k)` (points), `s(k-1)` (slopes `(f(i+1)-f(i))/(x(i+1)-x(i))`), `ierr`.

- **`FUNCTION pchst(arg1, arg2) RESULT(sign_product)`**:
  - **Description:** A utility function that determines if `arg1` and `arg2` have the same sign, opposite signs, or if one is zero. Returns `+1.0`, `-1.0`, or `0.0` respectively. It avoids direct multiplication to prevent overflow/underflow.
  - **Arguments:** `arg1, arg2 :: REAL(KIND=DP), INTENT(IN)`.

# Important Variables/Constants

The module's behavior is primarily controlled by the input arguments to its subroutines, especially:
- `x(:), f(:)`: The input data points.
- `d(:)`: The calculated or input derivative values at the data points.
- `spline :: LOGICAL` in `pchez`: Selects between PCHIP and cubic spline.
- `ic(:), vc(:)` in `pchsp`: Define boundary conditions for splines.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP` (double precision kind parameter).
- The module is largely self-contained for its mathematical logic.
- It is used by other modules that require interpolation of tabulated data. For example, `KEDF_DenDec` uses `pchez` and `pchev` to interpolate atomic core densities read from files.

**Typical Workflow:**
1.  Call `pchez` with the data `(x, f)` and a choice of `spline` to compute the derivative array `d`. This sets up the interpolant.
2.  Call `pchev` with `(x, f, d)` and new points `xval` to get interpolated function values `fval` and derivative values `dval` at these new points.
