# Overview

The `RhoOptN` module provides the main driver subroutine, `ProjNormRhoMinimization`, for optimizing the electron density variable `phi` (where `phi^2 = rho`, the electron density) while strictly conserving the total number of electrons. This electron number conservation is achieved by parameterizing trial steps using a mixing angle `theta`, such that `phi_new = phi_old * cos(theta) + dir_transformed * sin(theta)`, where `dir_transformed` is a search direction made orthogonal to `phi_old` and then normalized.

The module can employ several different algorithms to determine the raw search direction (`dir_type` parameter):
1.  Truncated Newton (TN): Solves `H*p = -g` using CG, where `H*p` is found by finite differences. (via `RhoDirNew`)
2.  Conjugate Gradient (CG): Uses Polak-Ribiere or Hager-Zhang formulas. (via `RhoDirCG`)
3.  Limited-memory BFGS (L-BFGS): Quasi-Newton method. (via `RhoDirBFGS`)

Once a raw direction is obtained, it's processed by routines from `RhoLineSearch` to get the initial `theta` and `dE/dtheta`, and then a line search for `theta` is performed (again, using routines from `RhoLineSearch`, which may involve `Dcsrch` or specialized methods for spin). The module also handles iterations for fixed magnetization calculations and special iteration loops for KEDFs like WGCD or EvW that require self-consistency for their parameters.

**References:**
- Press, W.H., et. al. Numerical Recipes.
- Gill, P.E., et.al. Practical Optimization.
- Jiang Hong and Weitao Yang, J. Chem. Phys. (2004) (for chemical potential formula).
- M.C. Payne, M.P. Teter, D.C. Allen, et al., Rev. Mod. Phys. (1992) (general DFT optimization context).

# Key Components

- **`MODULE RhoOptN`**: The main container module.
  - `infoUnit :: INTEGER`: File unit for detailed iteration output (default 6, screen).

- **`SUBROUTINE ProjNormRhoMinimization(energy, dir_type)`**:
  - **Description:** The primary driver for electron density optimization.
    1.  Initializes `phi = sqrt(rhoR)`.
    2.  Calculates initial `dE/dPhi` (using `CalculatePotentialPlus`) and `dL/dPhi = dE/dPhi - 2*mu*phi` (after calculating `mu` via `ChemicalPotential`).
    3.  Prints headers and initial state information.
    4.  Enters the main optimization loop (`DensityOptimization` internal subroutine).
    5.  Handles special iteration loops for WGCD (`kinetic == 12`) or EvW (`kloc > 0`) KEDFs, where `ComputeFrMatrix` or `aloc` updates might trigger further density optimization cycles (`GOTO 10`).
    6.  Handles iterations for fixed magnetization (`fixedmag` and `numSpin == 2`).
    7.  After convergence or reaching limits, updates `rhoR = phi^2`.
    8.  Prints final energy and timing information.
  - **Arguments:**
    - `energy(:) :: REAL(KIND=DP), INTENT(OUT)`: Final energy components.
    - `dir_type :: INTEGER, INTENT(INOUT)`: Specifies the search direction algorithm (1: TN, 2: CG, 3: L-BFGS). Can be modified internally (e.g., if TN fails).

- **Internal Subroutines to `ProjNormRhoMinimization`:**
    - **`SUBROUTINE DensityOptimization`**:
      - **Description:** Contains the main loop for iteratively optimizing `phi`.
        1.  Calls `getNextDirection` to get a raw search direction (`dirNext`) based on `dir_type`.
        2.  Calls `getNextGradient` (from `RhoLineSearch`) to make `dirNext` orthogonal to `phi`, normalize it, and get initial `theta` and `gradth` (`dE/dtheta`).
        3.  Calls `checkGradient` (from `RhoLineSearch`) to ensure `gradth` indicates a descent direction; may switch to steepest descent if not.
        4.  Calls `LineSearchTN` (from `RhoLineSearch`) to find the optimal `theta` that satisfies Wolfe conditions (or similar). This involves iterative calls to `CalculatePotentialPlus` and `gradOftheta`.
        5.  Updates `phi_new = phi_old*cos(theta) + dir_transformed*sin(theta)`.
        6.  Updates `dEdPhi` and `dLdPhi` using the new `phi` and recalculated `mu`.
        7.  Prints step information via `StepInformation`.
        8.  Checks for convergence via `CheckExit`.
        9.  If L-BFGS (`dir_type==3`), calls `UpdateBFGS`.
        10. Updates `phi` and `dLdPhi_old`, `dirNext_old` for the next iteration.
    - **`SUBROUTINE getNextDirection`**:
      - **Description:** A simple dispatcher that calls the appropriate search direction routine based on `dir_type`:
        - `dir_type = 1`: Calls `NewtonDirection` (from `RhoDirNew`).
        - `dir_type = 2`: Calls `CGDirection` (from `RhoDirCG`).
        - `dir_type = 3`: Calls `BFGSDirection` (from `RhoDirBFGS`).
    - **`SUBROUTINE StepInformation(...)`**:
      - **Description:** Gathers various statistics (min/max density, min/max potential, electron number deviation) and calls `MinimizerReportSteps` to print them.
    - **`SUBROUTINE CheckExit(...)`**:
      - **Description:** Checks termination criteria for the `DensityOptimization` loop:
        - Potential norm `potNorm` vs. `pot_tol`.
        - Energy change `abs(E_new - E_old)` vs. `tole` over a few steps.
        - Maximum iterations `maxIter`.
        Sets `rhoOutcome` and `exit_flag`.
    - **`SUBROUTINE DenOptReport(dir_type)`**:
      - **Description:** Prints a specific header line to the output indicating which optimization method (TN, CG, BFGS) is being used.

# Important Variables/Constants

This module heavily relies on parameters imported from `RhoOptimizers`, including:
- `maxIter, lineMaxIter, pot_tol, tole, conv_check, fixedmag, niter_extra, nmag, tolm`.
And routines/parameters from `RhoLineSearch`, `RhoDirCG`, `RhoDirNew`, `RhoDirBFGS`.

Key local arrays:
- `phi, tempPhi`: Store `sqrt(rho)` for current and trial steps.
- `dirNext, dirNext_old`: Store current and previous search directions.
- `dEdPhi, dLdPhi, dLdPhi_old`: Store `dE/dPhi` and `dL/dPhi` (gradient of Lagrangian `L=E-mu*N`).
- `mu0(:)`: Chemical potential(s).
- `theta(:)`: Angle(s) for line search.
- `gradth(:)`: `dE/dtheta`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`RhoOptimizers`**: For many control parameters, flags, and counters.
- **`OutputFiles`, `Output`**: For `outputUnit`, `infoUnit`, `QUIT`, `WrtOut`, and reporting subroutines from `Report`.
- **`CellInfo`**: For `kinetic` flag (for WGCD/EvW special loops), `numSpin`, `cell%dV`, `cell%vol`.
- **`CONSTANTS`**: For `DP`, `bohr`.
- **`RhoLineSearch`**: For `gradOftheta`, `getNextGradient`, `checkGradient`, `LineSearchTN`.
- **`RhoDirCG`**: For `InnerProd` and `CGDirection` (if `dir_type=2`).
- **`RhoDirNew`**: For `NewtonDirection` (if `dir_type=1`).
- **`RhoDirBFGS`**: For `InitializeBFGS`, `CleanBFGS`, `UpdateBFGS`, `BFGSDirection`, `ChemicalPotential` (if `dir_type=3`).
- **`Timer`**: For `TimerStart`, `TimerStop`, `stopwatch`.
- **`MathFunctions`**: For `MinMaxVal`.
- **`Sys`**: For `rhoR` (initial density, updated at end) and `rho0` (used by WGCD/EvW).
- **`KEDF_EvW`**: For `kloc`, `aloc`, `tolk` (parameters for EvW self-consistency).
- **`KEDF_WGCD`**: For `InitializeWGCD`, `ComputeFrMatrix` (for WGCD self-consistency).
- **`CalPotPlus`**: `CalculatePotentialPlus` is called by `LineSearchTN` (via `RhoLineSearch`) to get `dE/dPhi` and energy at trial `phi` values.
- **`KEDF_HC10`**: For `nref` (used in `StepInformation` for output).
- **`KEDF_DenDec`**: For `do_den_dec` and `potDD` (if density decomposition is used, `potDD` is updated with `dLdPhi/(2*phi)`).

The `ProjNormRhoMinimization` subroutine is a high-level driver for various gradient-based density optimization methods. It ensures electron number conservation through its specific line search parameterization (`phi_new = phi_old*cos(theta) + ...`). It's the workhorse for finding the ground state electron density.
