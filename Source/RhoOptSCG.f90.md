# Overview

The `RhoOptSCG` module implements electron density optimization using either Steepest Descent (SD) or Conjugate Gradient (CG) methods, specifically applied to `phi = sqrt(rho)` to ensure density positivity and facilitate electron number conservation. The "SCG" in the module name likely refers to "Sqrt Conjugate Gradient."

The core routine `SqrtGradientMinimization` iteratively updates `phi`. In each iteration:
1.  A search direction (`dirNext`) is determined:
    - For Steepest Descent, `dirNext = -gradient`.
    - For Conjugate Gradient (Polak-Ribiere like), `dirNext = -gradient + gamma * previous_direction`.
2.  A line search (`BrentLineMin`) is performed along `dirNext` to find an optimal step size. Brent's method, a robust technique combining golden section search and parabolic interpolation, is used.
3.  The `phi` variable is updated using this step size, and the process repeats.

The line search and gradient calculations are performed with respect to a Lagrangian `L = E - mu * N_electrons` to conserve the total number of electrons. The `Step` function used within the line search also enforces this conservation by rescaling the trial density.

# Key Components

- **`MODULE RhoOptSCG`**: The main container module.
  - `sumRho :: REAL(KIND=DP)`: Module-level variable storing the integrated total electron density (`N_electrons`), used for normalization in the `Step` function.
  - `tiny2 :: REAL(KIND=DP)`: A small number, slightly less than `tiny` from `Constants`, used for numerical comparisons.

- **`SUBROUTINE SqrtGradientMinimization(grid, energy, cgmin, maxRhoIter)`**:
  - **Description:** The primary driver for Steepest Descent or Conjugate Gradient optimization of `phi = sqrt(rhoR)`.
    1.  Initializes `phi = SQRT(rhoR)` and calculates initial `sumRho`.
    2.  Calls `CalculatePotentialPlus` to get initial `dE/dphi` and `energy`.
    3.  Calls the internal `Constrain` function to get the initial constrained gradient `pot = dL/dphi`.
    4.  Enters the main optimization loop (up to `maxRhoIter` steps):
        a.  Checks for convergence based on `potentialNorm` vs. `pot_tol` or energy change vs. `tole`.
        b.  Determines the search direction `dirNext`:
            - If first iteration or `cgmin` is `.FALSE.` or a restart condition is met (line search failure): `dirNext = -pot` (Steepest Descent).
            - If `cgmin` is `.TRUE.`: Calculates Polak-Ribiere `gamma = MAX(0, SUM[pot*(pot-g)]/SUM[g*g])`, then `dirNext = -pot + gamma * dirNext_old`. (`g` is the previous `pot`).
        c.  Calls the internal `BrentLineMin` subroutine to perform a line search along `dirNext`, which finds an optimal step `lastStep` and updates `rhoR` (which is `phi`) and `finalEnergy` (which is `energy`).
        d.  Updates `g = pot` (stores old gradient).
        e.  Recalculates `dE/dphi` (via `CalculatePotentialPlus`) and then `pot = dL/dphi` (via `Constrain`) using the new `rhoR` (`phi`).
        f.  Reports step statistics.
    5.  After the loop, converts `rhoR` back from `phi` to actual density (`rhoR = rhoR**2`).
  - **Arguments:**
    - `grid :: TYPE(gridPack), INTENT(INOUT)`: Grid pack type (seems to primarily pass `rhoR` through `Sys` module).
    - `energy(:) :: REAL(KIND=DP), INTENT(OUT)`: Final energy components.
    - `cgmin :: LOGICAL, INTENT(IN)`: If `.TRUE.`, use Conjugate Gradient; otherwise, use Steepest Descent.
    - `maxRhoIter :: INTEGER, INTENT(IN)`: Maximum allowed iterations.

- **Internal Functions/Subroutines to `SqrtGradientMinimization`:**
    - **`FUNCTION Constrain(step, rho, mask, mu, dimX, dimY, dimZ, nspin) RESULT(constrained_step)`**:
      - **Description:** Calculates the Lagrange multiplier `mu = SUM(step*mask) / (2*SUM(rho*mask))` and returns `constrained_step = step - 2*rho*mu`. `step` is typically `dE/dphi` and `rho` is `phi`. This projects `dE/dphi` to get `dL/dphi`.
    - **`SUBROUTINE BrentLineMin(grid, unitDirection, lambda, finalEnergy, timeB, success)`**:
      - **Description:** Implements Brent's line minimization method to find an optimal step size `timeB` along `unitDirection`. It first calls `MinBracket` to bracket a minimum, then uses a combination of golden section search and parabolic interpolation to refine it. The objective function is `E(phi_trial)`, where `phi_trial` is constructed using the internal `Step` function.
      - **Arguments:** `grid` (for `rhoR`), `unitDirection` (`dirNext`), `lambda` (chemical potential from `Constrain`), `finalEnergy` (INOUT), `timeB` (OUT, optimal step), `success` (OUT).
    - **`SUBROUTINE MinBracket(...)`**:
      - **Description:** Attempts to bracket a minimum for `BrentLineMin`. Starts with points `(0, E_current)` and `(init_step, E_trial)`, then expands or contracts the interval to find `a,b,c` such that `E(b) < E(a)` and `E(b) < E(c)`.
    - **`FUNCTION Step(sqrtRho, time, dir, mask) RESULT(rho_new_step)`**:
      - **Description:** Calculates `phi_trial = sqrtRho_old + time*dir`, then `rho_trial = phi_trial**2`. It then rescales `rho_trial` to ensure `SUM(rho_trial * dV)` equals the module-level `sumRho` (total number of electrons).
      - **Return Value:** `rho_new_step` (the new, norm-conserved density `rho_trial`).

# Important Variables/Constants

- **Module-Level:**
    - `sumRho :: REAL(KIND=DP)`: Stores the total number of electrons.
    - `tiny2 :: REAL(KIND=DP)`: A small constant.
- **Imported from `RhoOptimizers`:** `maxIter` (used as `maxRhoIter`), `lineMaxIter`, `pot_tol`, `tole`, `rhoOutcome`, `lineTolPercent`, `cheatPot`.
- **Key Local Arrays in `SqrtGradientMinimization`:** `g` (previous gradient), `dirNext`, `pot` (current constrained gradient).
- **Line Search (`BrentLineMin`):** `timeA, timeB, timeC` (bracketing points), `energyA, energyB, energyC` (energies at these points).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `tiny`, `DP`.
- **`RhoOptimizers`**: For various control parameters and status variables.
- **`KEDF_WGCkernel`**: For `firstOrderWGC` flag (used in convergence check logic).
- **`OUTPUT`, `OutputFiles`**: For `outputUnit`, `WrtOut`.
- **`CellInfo`**: For `cell%dV` (used in `InnerProd` via `Constrain`'s calculation of `mu`, which is then used by `CalculatePotentialPlus` context, not directly in `Constrain`).
- **`MathFunctions`**: For `volume` (not directly used, but `cell%vol` from `CellInfo` is related), `MinMaxVal`.
- **`Sys`**: For `rhoR` (the primary density array being optimized, by operating on its square root) and `interior` mask. `gridPack` type is an argument to `SqrtGradientMinimization` but seems unused.
- **`Timer`**: For `TimerStart`, `TimerStop`, `stopwatch`.
- **`Report`**: For `MinimizerReportHeader`, `MinimizerReportSteps`, `MinimizerReportFooter`.
- **`CalPotPlus`**: `CalculatePotentialPlus` is called by `BrentLineMin` (and at the start of `SqrtGradientMinimization`) to get `dE/dphi` and energy at trial `phi` values.
- **`MPI_Functions`**: For `ReduceRealLevel1`, `mpiErr`, `rankGlobal`.

The `SqrtGradientMinimization` routine is a density optimizer. It uses either Steepest Descent or CG for the direction, and Brent's method (via `BrentLineMin`, `MinBracket`, and `Step`) for a norm-conserving line search on `phi = sqrt(rho)`. The `Constrain` function is used to calculate the effective gradient `dL/dphi`.
