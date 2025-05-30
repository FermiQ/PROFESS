# Overview

The `LnSrch_MOD` module provides a custom line search subroutine named `LnSrch`. This routine is designed to be used within iterative optimization algorithms, such as Conjugate Gradient or Quasi-Newton methods, to find an optimal step length `alpha` along a given search direction. The goal is to find an `alpha` that satisfies the strong Wolfe conditions, ensuring both sufficient decrease in the objective function and a sufficient decrease in the magnitude of its gradient.

The `LnSrch` subroutine appears to use a strategy based on bracketing the minimum and then using linear interpolation of the gradient `g = df/d(alpha)` (where `f` is the objective function) to predict a new trial `alpha` where `g` would be zero. It manages a set of points (`points(2,3)`) that define the current search interval or best guesses.

# Key Components

- **`MODULE LnSrch_MOD`**: The main container module.
  - **Module-level variables:**
    - `dsave(13) :: REAL(KIND=DP)`: An array used to save values from previous steps in the line search. `dsave(1)` stores the function value `f` at `alpha=0` (i.e., `f_initial`), and `dsave(2)` stores the gradient `g` at `alpha=0` (i.e., `g_initial`). The other elements are not explicitly used in the provided code. This array's name and size are reminiscent of the `SAVE` array in MINPACK's `DCSRCH`.
    - `points(2,3) :: REAL(KIND=DP)`: Stores information about two key points in the line search, `(alpha, f_at_alpha, g_at_alpha)`. `points(1,:)` typically represents the "left" or starting point of an interval, and `points(2,:)` represents the "right" or current best trial point.

- **`SUBROUTINE LnSrch(alp, f, g, task, c1, c2, xtol, stpMin, stpMax)`**:
  - **Description:** The core line search routine. It takes the current step length `alp`, function value `f`, and gradient `g` (derivative `df/dalp`) as input. The `task` string variable controls the state and output of the routine.
    - **Initialization (`task == "START"`):** Stores initial `f` and `g` in `dsave` and `points(1,:)`. Sets an initial trial `alp`. Sets `task = "FG"`.
    - **Iteration (`task == "FG"` or `"WARNING..."`):**
      1.  Checks Wolfe conditions:
          - Sufficient decrease: `f <= f_initial + c1*alp*g_initial`
          - Strong curvature: `abs(g) <= c2*abs(g_initial)`
          If both are met, sets `task = "CONV"` and returns.
      2.  If not converged, it updates the `points` array based on the current `alp`, `f`, and `g` to maintain a bracket or refine the search region.
      3.  Calls the internal function `newAlp` to propose a new step length `alp` using linear interpolation of gradients, aiming for where the interpolated gradient is zero.
      4.  Handles various scenarios like interpolation errors or `alp` going out of bounds `[stpMin, stpMax]`, potentially adjusting `alp` heuristically (e.g., `alp = alp * 5.0` or bisection) and setting `task` to a "WARNING" or "ERROR" state.
      5.  If a new valid trial `alp` is found, sets `task = "FG"` to request function/gradient evaluation at this new `alp`.
  - **Arguments:**
    - `alp :: REAL(KIND=DP), INTENT(INOUT)`: Current/New trial step length.
    - `f :: REAL(KIND=DP), INTENT(IN)`: Function value at `alp`.
    - `g :: REAL(KIND=DP), INTENT(IN)`: Gradient `df/dalp` at `alp`.
    - `task :: CHARACTER(LEN=70), INTENT(INOUT)`: Communication string indicating state or requesting action.
    - `c1, c2 :: REAL(KIND=DP), INTENT(IN)`: Tolerance parameters for the Wolfe conditions (typically `c1` is small, e.g., 1e-4, and `c2` is larger, e.g., 0.1-0.9).
    - `xtol :: REAL(KIND=DP), INTENT(IN)`: Tolerance for the width of the interval in `alpha` (a comment notes this is not currently respected).
    - `stpMin, stpMax :: REAL(KIND=DP), INTENT(IN)`: Minimum and maximum allowed values for `alp`.

- **`FUNCTION newAlp(x1, y1, x2, y2, stpMin, stpMax) RESULT(alpha_trial)`** (Internal to `LnSrch`):
  - **Description:** Calculates a new trial step length `alpha_trial` by linear interpolation. Given two points on the `g(alpha)` curve, `(x1, y1)` and `(x2, y2)` (where `x` is `alpha` and `y` is `g`), it finds the `alpha` where the line connecting these points crosses zero: `alpha_trial = -b / slope`, where `slope = (y1-y2)/(x1-x2)` and `b = y1 - x1*slope`.
  - **Bounds Check:** If the calculated `alpha_trial` is outside `[stpMin, stpMax]`, it returns `-1.0` to signal an issue.
  - **Arguments:** `x1, y1, x2, y2, stpMin, stpMax` (all `REAL(KIND=DP)`).

# Important Variables/Constants

- **Wolfe Condition Parameters:** `c1`, `c2`.
- **Step Length Bounds:** `stpMin`, `stpMax`.
- **Saved State:** `dsave` (initial f and g), `points` (current bracketing/trial points `(alpha, f, g)`).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP` (double precision kind parameter).
- **`MPI_Functions`**: For `rankGlobal`, used to restrict print statements to rank 0 in parallel execution.

**Workflow:**
An optimization routine (e.g., `IonOptCG2`) would typically use `LnSrch` as follows:
1.  Initialize: `task = "START"`, provide an initial guess for `alp`.
2.  Loop:
    a.  Call `LnSrch(alp, f_current, g_current, task, ...)`.
    b.  If `task == "FG"`, the optimizer evaluates the objective function and its gradient at the new `alp` returned by `LnSrch`. Update `f_current` and `g_current`. Loop back to (a).
    c.  If `task == "CONV"`, the line search has succeeded. Exit loop.
    d.  If `task == "ERROR..."` or other warning/failure, handle appropriately (e.g., stop optimization or try a different strategy).

The `LnSrch` routine itself manages the process of refining `alp` by updating its internal state (`points` array) and using the `newAlp` function to propose candidates.
