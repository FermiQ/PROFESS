# Overview

The `RhoOptSTN` module implements an electron density optimization algorithm based on the Truncated Newton method, specifically applied to `phi = sqrt(rho)` to intrinsically maintain positivity of the density `rho`. The "STN" likely stands for Sqrt Truncated Newton.

The core idea of a Truncated Newton method is to approximately solve the Newton equation `H * p = -g` (where `H` is the Hessian matrix, `p` is the search direction, and `g` is the gradient) using an iterative method like Conjugate Gradient (CG). The CG iteration is "truncated" before full convergence, yielding an approximate Newton direction. A line search is then performed along this direction.

This module includes:
- The main driver `SqrtNewtonMinimization`.
- An internal CG solver (`NewtonDirection`) where the Hessian-vector product `H*p` is calculated using finite differences.
- A line search mechanism using the external `Dcsrch` routine.
- Constraint handling (`Constrain`) to ensure operations are performed with respect to the Lagrangian `L = E - mu * N_electrons`, effectively making the search direction norm-conserving (electron number conserving).
- A `Step` function to update the density while preserving the total electron number.

# Key Components

- **`MODULE RhoOptSTN`**: The main container module.
  - `sumRho :: REAL(KIND=DP)`: A module-level variable storing the total number of electrons (integral of density), used for normalization in the `Step` function.

- **`SUBROUTINE SqrtNewtonMinimization(energy)`**:
  - **Description:** The primary driver for the Sqrt Truncated Newton optimization.
    1.  Initializes `phi = SQRT(rhoR)`. Calculates initial `sumRho`.
    2.  Calls `CalculatePotentialPlus` to get initial `dE/dphi` and `energy`.
    3.  Calls `Constrain` to get the initial constrained potential `pot = dL/dphi = dE/dphi - mu*phi`.
    4.  Enters the main optimization loop:
        a.  Checks for energy convergence (`tole`) or potential norm convergence (`pot_tol`).
        b.  Handles `cheatPot` logic (to initially use simpler KEDFs if `firstOrderWGC` allows).
        c.  If issues arise (line search failure and already in steepest descent, or stagnation), may switch to steepest descent (`dirNext = -pot`) or exit.
        d.  Otherwise, calls the internal `NewtonDirection` function to get a search direction `dirNext` by iteratively solving `H * p = -pot` using CG.
        e.  Performs a line search using `Dcsrch` along `dirNext`. The function being minimized is `E(phi_trial)`, and the gradient projected along `dirNext` is `SUM(pot_trial * dirNext)`.
            i.  Inside the `Dcsrch` loop, the internal `Step` function calculates `rho_trial = (phi_old + stp*dirNext)^2`, ensuring `SUM(rho_trial)` remains `sumRho`.
            ii. `CalculatePotentialPlus` gets `dE/dphi_trial` and `energy_trial` for `rho_trial`.
            iii. `Constrain` calculates `pot_trial = dL/dphi_trial`.
        f.  Updates `rhoR` (actually `phi`) based on the successful step `stp` from `Dcsrch`.
        g.  Updates `pot` with the new gradient.
        h.  Reports step statistics.
    5.  After the loop, updates `rhoR = rhoR**2` (as `rhoR` was storing `phi`).
  - **Arguments:** `energy(:) :: REAL(KIND=DP), INTENT(OUT)`: Final energy components.

- **Internal Functions/Subroutines to `SqrtNewtonMinimization`:**
    - **`FUNCTION Constrain(step, sqrtRho, mask, mu, dimX, dimY, dimZ, nspin) RESULT(constrained_step)`**:
      - **Description:** Calculates the Lagrange multiplier `mu = SUM(step) / SUM(sqrtRho)` (sum over active grid points defined by `mask`) and returns `constrained_step = step - sqrtRho * mu`. This projects `step` (representing `dE/dphi`) onto a space that is effectively orthogonal to `sqrtRho`, resulting in `dL/dphi`.
      - **Arguments:** `step` (`dE/dphi`), `sqrtRho` (`phi`), `mask` (`interior`), `mu` (OUT).
    - **`FUNCTION NewtonDirection(b, b2, lambda, iter) RESULT(x)`**:
      - **Description:** Solves the linear system `A*x = -b` using a Conjugate Gradient (CG) method, where `x` is the Newton direction (`dirNext`), `b` is the current constrained gradient (`dL/dphi` from `Constrain`), and `A` is the Hessian. The Hessian-vector product `A*p` (where `p` is the CG search direction) is approximated by a finite difference:
        `Ap = ( (dL/dPhi)_trial - b2 ) / epsilon`, where `(dL/dPhi)_trial` is evaluated at `phi + epsilon*p`, and `b2` is `dL/dPhi` at `phi`. `lambda` here is the `mu` from `Constrain`.
      - **Arguments:** `b` (`dL/dphi`), `b2` (also `dL/dphi`), `lambda` (`mu`), `iter` (OUT, CG iterations).
      - **Return Value:** `x` (the Newton direction).
    - **`FUNCTION Step(sqrtRho, time, dir, mask) RESULT(rho_new_step)`**:
      - **Description:** Takes a step `time` (alpha from line search) along direction `dir` from `sqrtRho_old`. Calculates `phi_trial = sqrtRho_old + time * dir`. Then, it computes `rho_trial = phi_trial**2`. Finally, it rescales `rho_trial` so that `SUM(rho_trial * dV)` equals the initial `sumRho` (total number of electrons), effectively conserving electron number.
      - **Arguments:** `sqrtRho` (`phi_old`), `time` (step length), `dir` (search direction), `mask` (`interior`).
      - **Return Value:** `rho_new_step` (the new density `rho_trial` after stepping and rescaling).

# Important Variables/Constants

- **Module-Level:**
    - `sumRho :: REAL(KIND=DP)`: Stores the total number of electrons, used for normalization in the `Step` function.
- **Imported from `RhoOptimizers`:**
    - `maxIter, lineMaxIter, pot_tol, tole, cheatPot, rhoOutcome`, etc.: Control parameters for the optimization loops.
- **Local to `SqrtNewtonMinimization`:**
    - `phi`: Stores `sqrt(rhoR)`.
    - `dirNext`: The search direction.
    - `pot`: The constrained gradient `dL/dphi`.
    - `lambda`: Lagrange multiplier / chemical potential calculated by `Constrain`.
    - `stp`: Step length found by `Dcsrch`.
- **External Routine:**
    - `Dcsrch`: A line search algorithm (e.g., from MINPACK) used to find an optimal step length `stp`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`RhoOptimizers`**: For various control parameters (tolerances, iteration limits, flags) and status variables.
- **`OUTPUT`, `OutputFiles`**: For `outputUnit`, `WrtOut`, and reporting subroutines from `Report`.
- **`CellInfo`**: For `cell%dV` (used in `InnerProd` via `NewtonDirection` and `CalculatePotentialPlus`).
- **`MPI_Functions`**: For `ReduceRealLevel1` (for global sums in `Constrain`, `NewtonDirection`, `Step`) and `rankGlobal`.
- **`Timer`**: For `TimerStart`, `TimerStop`, `stopwatch`.
- **`Report`**: For `MinimizerReportHeader`, `MinimizerReportSteps`, `MinimizerReportFooter`.
- **`KEDF_WGCkernel`**: For `firstOrderWGC` flag (used for `cheatPot` logic).
- **`Sys`**: For the initial `rhoR` and the `interior` mask.
- **`CalPotPlus`**: For `CalculatePotentialPlus`, used to get `dE/dphi` and energy at trial densities during line search and finite difference for Hessian-vector product.
- **`Hartree`**: For `SqrtPreconditioner` (its usage in `NewtonDirection` is currently bypassed).
- **`CONSTANTS`**: For `DP`, `machPrec`.
- **`Dcsrch` (External)**: The line search routine used within `SqrtNewtonMinimization`.

The `SqrtNewtonMinimization` routine is a sophisticated density optimizer. It combines a Truncated Newton approach (using an inner CG to solve Newton's equations) with a line search and careful handling of constraints (electron number conservation via the `Constrain` and `Step` functions working with `sqrt(rho)`). The finite difference approximation of Hessian-vector products makes it matrix-free.
