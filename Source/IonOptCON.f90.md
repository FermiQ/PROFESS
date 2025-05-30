# Overview

The `IonOptCG` module (the filename `IonOptCON.f90` might be a slight misnomer or an older naming convention, as the content implements a Conjugate Gradient method) provides procedures for optimizing the positions of ions within a simulation cell. The goal is to find the ionic configuration that minimizes the total energy of the system, which corresponds to a state where forces on the ions are minimized.

This module implements a Conjugate Gradient (CG) algorithm. The specific variant of CG seems to be Polak-Ribiere-like, based on the calculation of the `gamma` parameter. A crucial part of this CG implementation is the line search, for which it utilizes the external `DCSRCH` routine (commonly found in MINPACK and other optimization libraries, designed by Moré and Thuente).

# Key Components

- **`MODULE IonOptCG`**: The main container for the CG ionic optimizer.

- **`SUBROUTINE GradientOptimization(Optimizer, rho, energy, forces, frozenIon)`**:
  - **Description:** This is the core routine that performs ionic relaxation using a Conjugate Gradient method.
    1.  It starts by calculating initial forces.
    2.  If not converged, it enters a loop:
        a.  Calculates the CG `gamma` parameter (Polak-Ribiere like: `MAX(SUM(F_new*(F_new-F_old))/SUM(F_old*F_old), 0.0)` where F is force).
        b.  Determines the new search direction `dirNext` using current forces and `gamma * previous_direction`.
        c.  Saves the current ion coordinates (`origCoord`).
        d.  Enters a line search loop using the `DCSRCH` subroutine. `DCSRCH` iteratively proposes a step length `stp`.
        e.  For each trial `stp`, ion positions are updated: `new_coords = MODULO(origCoord + cartesian_step, 1.0)`.
        f.  `RefreshIonTerms` is called, then `Optimizer` (the external electron density optimizer) is called, followed by `CalculateForces`.
        g.  The total energy `energy(1)` and the projection of forces along `dirNext` are fed back to `DCSRCH`.
        h.  The line search continues until `DCSRCH` signals convergence or an error.
        i.  After the line search, convergence of forces is checked.
        j.  The loop continues until overall force convergence or `maxIonStep` is reached.
  - **Arguments:**
    - `Optimizer :: EXTERNAL`: An external subroutine that optimizes the electron density for the current ion positions.
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: The electron density in real space.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: An array of energy components; `energy(1)` is used as the total energy for the line search.
    - `forces :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output array for calculated forces. `forces(:,:,1)` is used as the total force.
    - `frozenIon :: LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN)`: A logical mask indicating which ionic degrees of freedom are fixed.

# Important Variables/Constants

- **CG Algorithm Control:**
    - `maxIonStep :: INTEGER`: (Imported from `IonOptimizers`) Maximum number of CG iterations.
    - `forceCutoff :: REAL(KIND=DP)`: (Imported from `IonOptimizers`) Convergence criterion based on the maximum force on any ion.
- **Line Search Parameters (for `DCSRCH`):**
    - `ftol :: REAL(KIND=DP)`: Tolerance for the sufficient decrease condition (Armijo). Default: 0.01.
    - `gtol :: REAL(KIND=DP)`: Tolerance for the curvature condition (Wolfe). Default: 0.1.
    - `xtol :: REAL(KIND=DP)`: Relative width of the interval of uncertainty. Default: 1.0E-10.
    - `stpmin, stpmax :: REAL(KIND=DP)`: Minimum and maximum allowed step lengths for the line search. Defaults: 1.0E-100, 1000.0.
    - `task :: CHARACTER(len=60)`: Communication string with `DCSRCH`.
    - `isave :: INTEGER, DIMENSION(2)`, `dsave :: REAL(KIND=DP), DIMENSION(13)`: Integer and real save arrays for `DCSRCH`.
- **CG State Variables:**
    - `g :: REAL(KIND=DP), DIMENSION(SIZE(forces,1), SIZE(forces,2))`: Stores the current forces (gradient).
    - `dirNext :: REAL(KIND=DP), DIMENSION(SIZE(forces,1), SIZE(forces,2))`: Current CG search direction.
    - `gamma :: REAL(KIND=DP)`: Parameter for updating the search direction.
    - `stp :: REAL(KIND=DP)`: Step length determined by `DCSRCH`.
    - `fp :: REAL(KIND=DP)`: Projection of the force along the search direction (`-SUM(forces(:,:,1)*dirNext)`).
- **External Line Search Routine:**
    - `DCSRCH :: EXTERNAL`: The line search subroutine (e.g., from MINPACK by Moré and Thuente).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision).
- **`OUTPUT`, `OutputFiles`**: For writing messages and standard output (`WrtOut`, `outputUnit`).
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`).
- **`MPI_Functions`**: Module is used, but no direct MPI calls are visible in this specific file. It might be that `WrtOut` or other utilities are MPI-aware through this module.
- **`CellInfo`**: For `cell` derived type (to access `cell%ionTable`, `cell%cellReal`).
- **`MathFunctions`**: For `Inverse` (to convert Cartesian steps to fractional).
- **`IonOptimizers`**: For configuration parameters (`maxIonStep`, `forceCutoff`) and timer variables (`watch`, `watch2`).
- **`Report`**: For formatted reporting of geometry optimization progress (`GeometryMinimizerReportHeader`, `GeometryMinimizerReportFooter`, `GeometryMinimizerReportSteps`).
- **`RefreshIons`**: For `RefreshIonTerms` (called after ion positions are updated).
- **`CalForces`**: For `CalculateForces` (to compute forces on ions at each trial step in the line search).

**External Library Dependency:**
- The module relies on an external line search routine named `DCSRCH`. This is a standard routine from libraries like MINPACK.

The `GradientOptimization` subroutine drives the ion relaxation. It calculates an initial set of forces. Then, in a loop, it computes a CG search direction, and calls `DCSRCH` to find an appropriate step length along this direction. `DCSRCH` itself will trigger recalculations of energy and forces (via `Optimizer` and `CalculateForces`) at trial points. Once `DCSRCH` finds a suitable step, ion positions are updated, and the overall convergence is checked.
