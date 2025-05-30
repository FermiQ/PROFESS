# Overview

The `IonOptBFGS` module implements ionic position optimization using the limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm. L-BFGS is a quasi-Newton method that approximates the inverse Hessian matrix to guide the search for a minimum energy configuration. It is well-suited for large-scale optimization problems as it avoids storing the full Hessian.

This module takes the electronic density, energy, and an external routine for optimizing the electron density (`RhoOptimizer`) as inputs. It iteratively updates ion positions to minimize forces, calling `RhoOptimizer` and force calculation routines at each step.

# Key Components

- **`MODULE IonOptBFGS`**: The main container for the L-BFGS ionic optimizer.

- **`SUBROUTINE IonLBFGS(RhoOptimizer, rho, energy, forces, frozenIon)`**:
  - **Description:** Performs ion relaxation using the L-BFGS algorithm. It initializes the L-BFGS parameters, converts ion coordinates to a 1D array `X` suitable for the L-BFGS routine, and then enters a loop where it calls the L-BFGS solver (`LB2`). The L-BFGS solver, in turn, calls back to `ComputeEnergyAndGradient` to get the current total energy (`F`) and the negative of the forces (`G`). The loop continues until force convergence or other L-BFGS termination criteria are met.
  - **Arguments:**
    - `RhoOptimizer :: EXTERNAL`: An external subroutine that optimizes the electron density for the current ion positions.
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: The electron charge density in real space.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: An array of energy components; `energy(1)` is taken as the total energy `F`.
    - `forces :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output array where the calculated forces on ions are stored. `forces(:,:,1)` is used as the total force.
    - `frozenIon :: LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN)`: A logical mask indicating which ionic degrees of freedom are frozen (not optimized).

  - **Internal Subroutine:**
    - **`SUBROUTINE ComputeEnergyAndGradient(relaxDone)`**:
      - **Description:** This subroutine is called by the L-BFGS algorithm (indirectly, as L-BFGS typically requires function and gradient evaluations). It first calls `RefreshIonTerms` to update ion-dependent quantities, then calls the `RhoOptimizer` to relax the electron density, and finally calls `CalculateForces` to get the forces on ions. The total energy `F` (from `energy(1)`) and the gradient `G` (negative of `forces(:,:,1)`) are updated. It also checks if the maximum force is below `forceCutoff` to set `relaxDone`.
      - **Arguments:**
        - `relaxDone :: LOGICAL, INTENT(OUT)`: Set to `.TRUE.` if force convergence is achieved.

# Important Variables/Constants

- **L-BFGS Algorithm Parameters:**
    - `N :: INTEGER`: Number of variables to optimize (3 * number of ions).
    - `M :: INTEGER`: Number of past corrections (gradient and position differences) stored by L-BFGS (default: 5).
    - `IPRINT(2) :: INTEGER`: Controls output from the L-BFGS routine (default: no output).
    - `DIAGCO :: LOGICAL`: Flag indicating whether a user-provided diagonal Hessian approximation is used (default: `.FALSE.`).
    - `EPS :: REAL(KIND=DP)`: Convergence tolerance for the gradient norm (forces) within L-BFGS itself (though PROFESS often uses its own `forceCutoff`).
    - `XTOL :: REAL(KIND=DP)`: Tolerance for step length in line search (default: 1.0D-4).
    - `W(:) :: REAL(KIND=DP), ALLOCATABLE`: Working array required by the L-BFGS algorithm.
    - `IFLAG :: INTEGER`: Communication flag between the caller and the L-BFGS routine.
- **Optimization Variables:**
    - `X :: REAL(KIND=DP), DIMENSION(3*cell%numIon)`: A 1D array holding the Cartesian coordinates of all ions, used by the L-BFGS routine.
    - `F :: REAL(KIND=DP)`: Current total energy of the system, passed to L-BFGS.
    - `G :: REAL(KIND=DP), DIMENSION(3*cell%numIon)`: The negative of the current forces on ions (gradient of energy), passed to L-BFGS.
- **Convergence Control:**
    - `forceCutoff :: REAL(KIND=DP)`: (Imported from `IonOptimizers`) The primary convergence criterion based on the maximum force on any ion.
    - `nMax :: INTEGER`: Maximum number of allowed calls to `ComputeEnergyAndGradient` (default: 1000).
- **External L-BFGS Routine:**
    - `LB2 :: EXTERNAL`: The name of the external L-BFGS subroutine that this module calls. This is not a Fortran module but an external library routine.
    - `/LB3/ :: COMMON MP,LP,GTOL,STPMIN,STPMAX`: A common block used to pass additional parameters to the `LB2` routine.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision).
- **`CellInfo`**: For `cell` derived type (to access `cell%numIon`, `cell%cellReal`, `cell%ionTable`).
- **`OUTPUT`**: For `WrtOut` (writing messages).
- **`OutputFiles`**: For `outputUnit` (standard output unit).
- **`REPORT`**: For formatted reporting of geometry optimization progress (`GeometryMinimizerReportHeader`, `GeometryMinimizerReportFooter`, `GeometryMinimizerReportSteps`).
- **`IonOptimizers`**: For `forceCutoff` (convergence threshold) and timer variables (`watch`, `watch2`).
- **`TIMER`**: For `TimerStart` and `TimerStop` (performance timing).
- **`MathFunctions`**: For `Inverse` (to convert fractional to Cartesian coordinates).
- **`MPI_Functions`**: For `message` (MPI-aware message writing).
- **`RefreshIons`**: For `RefreshIonTerms` (called after ion positions change).
- **`CalForces`**: For `CalculateForces` (to compute forces on ions).

**External Library Dependency:**
- The module critically depends on an external L-BFGS library routine, referred to as `LB2`. The communication with this routine also involves a Fortran `COMMON` block (`/LB3/`). This is a significant dependency that is not managed through standard Fortran module `USE` statements.

The `IonLBFGS` subroutine orchestrates the optimization. It prepares data in the format required by the `LB2` routine and then enters a loop. Inside the loop, `LB2` proposes new ion positions. `IonLBFGS` then updates the actual ion coordinates, calls `ComputeEnergyAndGradient` (which involves electron density relaxation via `RhoOptimizer` and force calculation), and feeds the new energy and forces back to `LB2`. This process repeats until convergence.
