# Overview

The `IonOptCG2` module provides subroutines for optimizing the positions of ions in a simulation cell to find a minimum energy configuration. It employs Conjugate Gradient (CG) methods, which are iterative optimization algorithms well-suited for multidimensional problems. The module specifically mentions support for Polak-Ribiere, Hager-Zhang, and Dai-Yuan variants of the CG algorithm.

A key part of the CG method is a line search procedure, and this module utilizes a custom line search routine (`lnsrch` from `LnSrch_MOD`) that aims to satisfy strong Wolfe conditions to ensure robust convergence.

# Key Components

- **`MODULE IonOptCG2`**: The main container for the CG ionic optimizer.

- **`SUBROUTINE ConjugateGradient(RhoOptimizer, rho, energy, forces, frozenIon)`**:
  - **Description:** This is the primary routine for performing ionic relaxation using a Conjugate Gradient method. It iteratively adjusts ion positions based on calculated forces. In each CG step, a search direction (`dir`) is determined using the current forces (gradient) and information from the previous step (previous gradient and search direction, combined via the `gam` parameter specific to the chosen CG variant). A line search (`lnsrch`) is then performed along this direction to find an optimal step length (`alpha`). Ion positions are updated, and the process repeats until forces converge below `forceCutoff` or the maximum number of steps (`maxIonStep`) is reached.
  - **Arguments:**
    - `RhoOptimizer :: EXTERNAL`: An external subroutine responsible for optimizing the electron density for the current set of ion positions.
    - `rho :: REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN)`: The electron density in real space.
    - `energy :: REAL(kind=DP), DIMENSION(:), INTENT(IN)`: An array of energy components; `energy(1)` is used as the total energy for line search and reporting.
    - `forces :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output array where calculated forces on ions are stored. `forces(:,:,1)` is used as the total force for CG.
    - `frozenIon :: LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN)`: A logical mask indicating which ionic degrees of freedom are fixed (frozen) and not subject to optimization.

  - **Internal Helper Subroutines:**
    - **`SUBROUTINE CheckExit(maximumForce, relaxDone)`**:
      - **Description:** Checks if the maximum force component among all ions (`maximumForce`) has fallen below the `forceCutoff` threshold. Sets `relaxDone` accordingly.
      - **Arguments:**
        - `maximumForce :: REAL(KIND=DP), INTENT(OUT)`: Calculated maximum force.
        - `relaxDone :: LOGICAL, INTENT(OUT)`: Flag indicating convergence.
    - **`SUBROUTINE StepInfo(step, maxForce, stepSizeAB, energy)`**:
      - **Description:** Prints information about the current CG iteration, including step number, maximum force, step size (`alpha`), and total energy.
      - **Arguments:** `step`, `maxForce`, `stepSizeAB` (alpha), `energy`.
    - **`SUBROUTINE gradOfAlpha(forces, gradient, numIons, gradAlp)`**:
      - **Description:** Computes the derivative of the total energy with respect to the line search step length `alpha` (i.e., `dE/d(alpha)`). This is calculated as `-SUM(forces .dot. search_direction)`. This gradient is crucial for the line search algorithm.
      - **Arguments:**
        - `forces :: REAL(KIND=DP), INTENT(IN)`: Current forces on ions.
        - `gradient :: REAL(KIND=DP), INTENT(IN)`: Current search direction (`dir`).
        - `numIons :: INTEGER, INTENT(IN)`: Number of ions.
        - `gradAlp :: REAL(KIND=DP), INTENT(OUT)`: Calculated `dE/d(alpha)`.

# Important Variables/Constants

- **CG Algorithm Control:**
    - `cg_type :: INTEGER`: (Imported from `IonOptimizers`) Selects the CG variant (1: Polak-Ribiere, 2: Hager-Zhang, 3: Dai-Yuan).
    - `maxIonStep :: INTEGER`: (Imported from `IonOptimizers`) Maximum number of CG iterations allowed.
    - `forceCutoff :: REAL(KIND=DP)`: (Imported from `IonOptimizers`) Convergence criterion based on the maximum force.
- **Line Search Parameters (for `lnsrch`):**
    - `ftol :: REAL(KIND=DP)`: Tolerance for the sufficient decrease condition (Armijo condition).
    - `gtol :: REAL(KIND=DP)`: Tolerance for the curvature condition (Wolfe condition).
    - `xtol :: REAL(KIND=DP)`: Minimum relative step length.
    - `stpmin, stpmax :: REAL(KIND=DP)`: Minimum and maximum allowed step lengths for the line search. `stpmax` can be dynamically adjusted based on `nnDist`.
- **CG State Variables:**
    - `dir :: REAL(KIND=DP), DIMENSION(cell%numIon,3)`: Current search direction.
    - `oldGrad, newGrad :: REAL(KIND=DP), DIMENSION(cell%numIon,3)`: Forces from the previous and current steps, respectively.
    - `gam :: REAL(KIND=DP)`: Parameter used to update the search direction in CG (e.g., Polak-Ribiere `beta`).
    - `alpha :: REAL(KIND=DP)`: Step length determined by the line search.
- **`nnDist :: REAL(KIND=DP)`**: (Imported from `NearestDistance`) Minimum allowed interatomic distance; used to cap `stpmax` in line search.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision).
- **`Output`, `OutputFiles`**: For writing messages and standard output (`WrtOut`, `outputUnit`).
- **`MPI_Functions`**: For `message` (MPI-aware message writing).
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`).
- **`CellInfo`**: For `cell` derived type (to access `cell%numIon`, `cell%cellReal`, `cell%ionTable`).
- **`MathFunctions`**: For `Inverse` (to convert fractional to Cartesian coordinates and vice-versa).
- **`LnSrch_MOD`**: Critically depends on this module for the line search algorithm (`lnsrch`) and for monitoring line search progress (`work` aliased to `points`).
- **`Report`**: For formatted reporting of geometry optimization progress (`GeometryMinimizerReportHeader`, `GeometryMinimizerReportFooter`, `GeometryMinimizerReportSteps`).
- **`IonOptimizers`**: For configuration parameters like `cg_type`, `maxIonStep`, `forceCutoff`, and timer variables (`watch`, `watch2`).
- **`NearestDistance`**: For `CheckNearestDistanceAtoms` and `nnDist` to prevent atoms from getting too close during optimization.
- **`RefreshIons`**: For `RefreshIonTerms` (called after ion positions are updated).
- **`CalForces`**: For `CalculateForces` (to compute forces on ions at each step).

The `ConjugateGradient` subroutine is the main driver. It initializes CG parameters, then enters a loop. In each iteration, it calculates the `gam` factor based on the chosen `cg_type`, determines a new search `dir`ection, performs a line search using `lnsrch` (which itself will call `RhoOptimizer` and `CalculateForces` via `gradOfAlpha` and `CheckExit` within its trial steps), updates ion positions with the found `alpha`, and checks for convergence.
