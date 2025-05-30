# Overview

The `MathSplines` module provides a set of routines for performing cubic spline interpolation. This package allows users to:
1.  Determine the second derivatives of a cubic spline at a given set of data points (knots), subject to specified boundary conditions. This is done by `spline_cubic_set`.
2.  Evaluate the resulting cubic spline, along with its first and second derivatives, at any arbitrary point within or outside (via extrapolation) the range of the original data points. This is done by `spline_cubic_val`.

The cubic spline is a piecewise cubic polynomial that ensures continuity of the function itself, its first derivative, and its second derivative across all knot points. This results in a smooth interpolating curve.

The code is attributed to John Burkardt's collection of Fortran routines and is licensed under the GNU Lesser General Public License (LGPL). It is likely based on well-established numerical algorithms for spline interpolation, possibly from sources like Carl de Boor or Fred Fritsch.

# Key Components

- **`MODULE MathSplines`**: The main container for the spline routines.

- **`SUBROUTINE spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp)`**:
  - **Description:** Computes the array of second derivative values (`ypp`) of the cubic spline at each knot `t(i)`. This is achieved by solving a tridiagonal system of linear equations derived from the condition that the first derivative of the spline must be continuous across interior knots, combined with user-specified boundary conditions.
  - **Arguments:**
    - `n :: INTEGER, INTENT(IN)`: The number of data points (knots). Must be at least 2.
    - `t(n) :: REAL(KIND=DP), INTENT(IN)`: The array of knot locations (abscissas), which must be strictly increasing.
    - `y(n) :: REAL(KIND=DP), INTENT(IN)`: The array of data values (ordinates) at each knot `t(i)`.
    - `ibcbeg :: INTEGER, INTENT(IN)`: Flag for the boundary condition at the left endpoint (`t(1)`):
        - `0`: Quadratic spline over the first interval (implies `ypp(1) = ypp(2)`).
        - `1`: The first derivative `y'(t(1))` is specified by `ybcbeg`.
        - `2`: The second derivative `y''(t(1))` is specified by `ybcbeg`.
    - `ybcbeg :: REAL(KIND=DP), INTENT(IN)`: The value for the left boundary condition if `ibcbeg` is 1 or 2.
    - `ibcend :: INTEGER, INTENT(IN)`: Flag for the boundary condition at the right endpoint (`t(n)`), similar to `ibcbeg`.
    - `ybcend :: REAL(KIND=DP), INTENT(IN)`: The value for the right boundary condition if `ibcend` is 1 or 2.
    - `ypp(n) :: REAL(KIND=DP), INTENT(OUT)`: The calculated second derivative values at each knot `t(i)`.

- **`SUBROUTINE spline_cubic_val(n, t, y, ypp, tval, yval, ypval, yppval)`**:
  - **Description:** Evaluates the cubic spline (defined by knots `t`, data values `y`, and second derivatives `ypp` from `spline_cubic_set`) at a specified point `tval`. It returns the interpolated function value `yval`, its first derivative `ypval`, and its second derivative `yppval` at `tval`. If `tval` is outside the range `[t(1), t(n)]`, extrapolation is performed using the cubic polynomial of the nearest end interval.
  - **Arguments:**
    - `n :: INTEGER, INTENT(IN)`: Number of data points.
    - `t(n) :: REAL(KIND=DP), INTENT(IN)`: The knot values.
    - `y(n) :: REAL(KIND=DP), INTENT(IN)`: The data values at the knots.
    - `ypp(n) :: REAL(KIND=DP), INTENT(IN)`: The second derivatives at the knots (from `spline_cubic_set`).
    - `tval :: REAL(KIND=DP), INTENT(IN)`: The point at which to evaluate the spline.
    - `yval :: REAL(KIND=DP), INTENT(OUT)`: Interpolated function value at `tval`.
    - `ypval :: REAL(KIND=DP), INTENT(OUT)`: Interpolated first derivative value at `tval`.
    - `yppval :: REAL(KIND=DP), INTENT(OUT)`: Interpolated second derivative value at `tval`.

- **`SUBROUTINE rvec_bracket(n, x, xval, left, right)`**:
  - **Description:** An internal helper routine that searches a sorted array `x` (typically the knots `t`) to find the interval `[x(left), x(right)]` that contains or is nearest to the given value `xval`. It uses a bisection-like search for efficiency.
  - **Arguments:**
    - `n :: INTEGER, INTENT(IN)`: Length of the array `x`.
    - `x(n) :: REAL(KIND=DP), INTENT(IN)`: The sorted array.
    - `xval :: REAL(KIND=DP), INTENT(IN)`: The value to be bracketed.
    - `left, right :: INTEGER, INTENT(OUT)`: Indices such that `x(left) <= xval <= x(right)` (or `xval` is outside the bounds).

- **`SUBROUTINE s3_fs(a1, a2, a3, n, b, x)`**:
  - **Description:** An internal helper routine that factors and solves a tridiagonal linear system `A*x=b`. The tridiagonal matrix `A` is defined by its sub-diagonal `a1` (from index 2 to `n`), diagonal `a2` (1 to `n`), and super-diagonal `a3` (1 to `n-1`). This is used by `spline_cubic_set` to solve for the unknown second derivatives `ypp`.
  - **Arguments:** `a1`, `a2`, `a3`, `n`, `b` (right-hand side, overwritten), `x` (solution vector).

# Important Variables/Constants

The module primarily consists of subroutines. The key data passed between them are:
- `t(:)`: The knot locations (x-values).
- `y(:)`: The data values at the knots (y-values).
- `ypp(:)`: The second derivatives of the spline at the knots. This array is computed by `spline_cubic_set` and used as input by `spline_cubic_val`.
- Boundary condition flags (`ibcbeg`, `ibcend`) and values (`ybcbeg`, `ybcend`).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- The module is largely self-contained for its spline interpolation logic.
- It is utilized by other modules in the PROFESS codebase that require interpolation of one-dimensional tabulated data. For example:
    - `KEDF_HC10` and `KEDF_WGCkernel` use `spline_cubic_set` and `spline_cubic_val` for interpolating kinetic energy kernel functions.
    - `LocalPseudoPot` uses logic very similar to `spline_cubic_val` for interpolating pseudopotentials, and relies on pre-calculated second derivatives (`psp%potDD`) which would have been generated by a routine like `spline_cubic_set`.

**Typical Workflow:**
1.  Provide a set of data points `(t(i), y(i))`.
2.  Call `spline_cubic_set` with `t`, `y`, and desired boundary conditions to calculate the second derivatives `ypp` at each knot. This "sets up" the spline.
3.  Call `spline_cubic_val` with `t`, `y`, `ypp`, and a new point `tval` to obtain the interpolated value `yval` (and optionally `ypval`, `yppval`) at `tval`.
