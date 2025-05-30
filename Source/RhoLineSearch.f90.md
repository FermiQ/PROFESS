# Overview

The `RhoLineSearch` module provides subroutines dedicated to performing line searches as part of an electron density (or `phi = sqrt(rho)`) optimization process. A key feature of the line search strategy implemented here is the explicit conservation of the total number of electrons. This is achieved by parameterizing the trial `phi` as:

`phi_new(theta) = phi_old * cos(theta) + dir_transformed * sin(theta)`

where `dir_transformed` is a search direction that has been made orthogonal to `phi_old` and then normalized. The line search then becomes an optimization problem for the angle `theta`. The module includes routines to:
1.  Prepare the initial transformed search direction and initial `theta` and `dE/dtheta`.
2.  Calculate `dE/dtheta` for a given `theta`.
3.  Perform the line search for `theta`, potentially using external routines like `Dcsrch` or custom methods for multi-spin cases.

**References:**
The module comments cite general optimization texts:
- Press, W.H., et. al. Numerical Recipes.
- Gill, P.E., et.al. Practical Optimization.

# Key Components

- **`MODULE RhoLineSearch`**: The main container module.

- **Module-Level Parameter:**
    - `maxLine :: INTEGER`: Maximum number of iterations allowed within a line search procedure (default: 200).

- **`SUBROUTINE getNextGradient(totEleNum, phi, dEdPhi, theta, gradth, dirNext)`**:
  - **Description:** Prepares the search for the optimal `theta`.
    1.  Orthogonalizes the input search direction `dirNext` with respect to the current `phi` for each spin channel: `dirNext_orth = dirNext - phi * <phi|dirNext> / <phi|phi>`.
    2.  Normalizes `dirNext_orth` (scaled by `sqrt(totEleNum)`).
    3.  Calculates an initial guess for `theta`. For `numSpin=1`, `theta = min(theta_current_guess, |dirNext_orth| / |phi + dirNext_original|)`. For `numSpin=2`, `theta` is initialized to 0.
    4.  Computes the initial gradient `gradth = dE/dtheta` at this `theta` by calling `gradOftheta`.
  - **Arguments:**
    - `totEleNum(:) :: REAL(KIND=DP), INTENT(IN)`: Total electrons per spin channel.
    - `phi(:,:,:,:) :: REAL(KIND=DP), INTENT(IN)`: Current `phi`.
    - `dEdPhi(:,:,:,:) :: REAL(KIND=DP), INTENT(IN)`: Current `dE/dPhi`.
    - `theta(:) :: REAL(KIND=DP), INTENT(OUT)`: Initialized/updated `theta` values.
    - `gradth(:) :: REAL(KIND=DP), INTENT(OUT)`: Initialized `dE/dtheta` values.
    - `dirNext(:,:,:,:) :: REAL(KIND=DP), INTENT(OUT)`: The input search direction, modified to be the orthogonalized and normalized direction.

- **`FUNCTION gradOftheta(theta, dEdPhi, phi, dirNext) RESULT(dEdTheta_val)`**:
  - **Description:** Calculates the derivative of the total energy `E` with respect to the mixing angle `theta`. The formula is `dE/dtheta = SUM_grid { (dE/dPhi_new) * d(Phi_new)/dtheta * dV }`, where `d(Phi_new)/dtheta = -phi_old*sin(theta) + dir_transformed*cos(theta)`. `dE/dPhi_new` is the gradient of energy with respect to the *new* `phi` that corresponds to the current `theta`.
  - **Arguments:**
    - `theta :: REAL(KIND=DP), INTENT(IN)`: The current angle.
    - `dEdPhi(:,:,:) :: REAL(KIND=DP), INTENT(IN)`: `dE/dPhi` evaluated at `phi_new(theta)`.
    - `phi(:,:,:) :: REAL(KIND=DP), INTENT(IN)`: `phi_old`.
    - `dirNext(:,:,:) :: REAL(KIND=DP), INTENT(IN)`: The (orthogonalized and normalized) transformed direction.
  - **Return Value:** `dEdTheta_val :: REAL(KIND=DP)`.

- **`SUBROUTINE checkGradient(switchStp, switchStp2, gradth, grad_exit_flag, rhoOutcome, dirNext, dLdPhi)`**:
  - **Description:** Checks the sign of `gradth` (dE/dtheta). If `gradth >= 0` at the start of a line search (theta=0), it implies the current search direction `dirNext` is not a descent direction. This routine then attempts to switch `dirNext` to the steepest descent direction (`-dLdPhi`). If `gradth` is still non-negative even with steepest descent, it signals a potential problem or convergence (`rhoOutcome = 1`).
  - **Arguments:** Various flags and gradients to control and report the state.

- **`SUBROUTINE LineSearchTN(iter, theta, energy, gradth, ftol, gtol, xtol, success, rhoOutcome, tempPhi, phi, dirNext, dEdPhi)`**:
  - **Description:** The main line search driver.
    - For `numSpin == 1`: It uses the external `Dcsrch` routine (a common line search algorithm, e.g., from MINPACK) to find `theta(1)` that satisfies Wolfe conditions. The function to minimize is effectively `E(theta)`, and its gradient is `gradth(1)`. Inside the `Dcsrch` loop, trial `tempPhi` values are generated, `CalculatePotentialPlus` is called to get new `energy` and `dEdPhi` at `tempPhi`, and then `gradOftheta` is called to get the new `gradth` for `Dcsrch`.
    - For `numSpin == 2`: It calls `SpinThetaSearch`.
  - **Arguments:** Many control and state variables for the line search and overall optimization.

- **`SUBROUTINE SpinThetaSearch(...)`**:
  - **Description:** Performs a line search for the 2-component vector `(theta_up, theta_down)` in spin-polarized cases. It appears to use an internal CG-like method to optimize this `theta` vector. The "objective function" for this inner optimization is implicitly related to minimizing the norm of the `gradth` vector, and its "gradient" is `d(E)/d(alpha_theta) = gradth_up * thetadir_up + gradth_down * thetadir_down`.
  - **Contains Internal Functions:**
    - `FUNCTION thetaCGDirection(...)`: Calculates a CG-like search direction for the `theta` vector.
    - `FUNCTION ThetaInnerProd(...)`: Calculates the dot product of two 2-component `theta` vectors.

# Important Variables/Constants

- **`maxLine :: INTEGER`**: Module-level parameter for maximum line search iterations (default: 200).
- **Line Search Parameters (for `Dcsrch` or similar):** `ftol`, `gtol`, `xtol`, `stpMin`, `stpMax`.
- **State Variables for Line Search:** `theta(:)`, `gradth(:)`, `energy(:)`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`, `PI`.
- **`MPI_Functions`**: For `ReduceRealLevel1` (used in `gradOftheta`), `Title`, `StartClock`, `StopClock`, `rankGlobal`.
- **`RhoDirCG`**: For `InnerProd` (used in `getNextGradient` for orthogonalization and normalization).
- **`CellInfo`**: For `numSpin`, `cell%dV`.
- **`Output`**: For `WrtOut`, `QUIT`, `message`.
- **`RhoOptimizers`**: For parameters `numEnergyLineMin`, `maxIter` (used in `LineSearchTN` to limit internal `Dcsrch` loop).
- **`CalPotPlus`**: For `CalculatePotentialPlus`, which is called repeatedly within the `LineSearchTN` (specifically, inside the `Dcsrch` loop or `SpinThetaSearch`) to evaluate energy and `dE/dPhi` at trial `phi` values corresponding to trial `theta`s.
- **`Dcsrch` (External Routine)**: The `LineSearchTN` subroutine explicitly calls `Dcsrch` for single-spin line searches. This is a common external Fortran routine for line searches satisfying Wolfe conditions (e.g., from MINPACK).

This module is a critical part of density optimizers that require electron number conservation (which is most of them). The optimizers (like `RhoOptN`) would first determine a raw search direction (e.g., from `RhoDirCG` or `RhoDirBFGS`), then call `getNextGradient` to transform it and get initial `theta` and `dE/dtheta`, and finally call `LineSearchTN` to find the optimal `theta`. The new `phi` is then constructed using this `theta`.
