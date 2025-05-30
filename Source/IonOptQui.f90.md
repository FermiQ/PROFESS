# Overview

The `IonOptQui` module implements an ionic position optimization algorithm known as "QuickMin." This method is a projected Velocity Verlet-like algorithm designed to find the minimum energy configuration of ions. It simulates a type of molecular dynamics where velocities are adjusted based on forces, but with specific rules to guide the system towards an energy minimum rather than conserving energy for true dynamics.

The core idea of QuickMin is to move ions based on current forces and velocities. If a step leads to an increase in energy or if the direction of movement (velocity) opposes the current forces too much (i.e., `forces .dot. velocity < 0`), the step is considered an "overstep." In such cases, the time step is reduced, and velocities might be reset. Otherwise, the time step is typically increased to accelerate convergence, and velocities are updated.

# Key Components

- **`MODULE IonOptQui`**: The main container for the QuickMin optimization routine.

- **`SUBROUTINE QuickMinOptimization(Optimizer, rho, energy, forces, frozenIon)`**:
  - **Description:** This is the primary routine that performs ionic relaxation using the QuickMin algorithm. It initializes velocities to zero and iteratively updates ion positions.
    1.  Calculates initial forces and energy.
    2.  Enters a loop that continues until forces converge below `forceCutoff` or `maxIonStep` is reached.
    3.  In each step:
        a.  Updates ion positions using current `velocity`, `forces`, and `timeStep` (similar to Velocity Verlet: `pos_new = pos_old + v*dt + 0.5*F/m*dt^2`, assuming mass m=1).
        b.  Calls `RefreshIonTerms` and then the external `Optimizer` to relax the electron density for the new ion positions.
        c.  Calculates new forces and energy (`testEnergy`).
        d.  Checks for "overstep": if `SUM(forces .dot. velocity) < 0` or if `testEnergy(1) > oldEnergy(1)`.
            -   If an overstep occurs (especially on the first step or if energy increases):
                -   The `timeStep` is significantly reduced (e.g., divided by 5).
                -   Velocities might be reset to zero.
                -   The ion positions are reverted to `lastCoord`.
                -   The loop continues to the next iteration with the adjusted parameters.
            -   If only a slight overstep where energy still decreases but `forces .dot. velocity < 0`:
                -   Velocities are reset to zero, but the new position is kept.
        e.  If not an overstep:
            -   The `timeStep` is increased (e.g., multiplied by 2) to accelerate.
            -   Velocities are updated: `velocity = (forces .dot. velocity_old) * forces / |forces|^2 + forces * timeStep` (projected update).
        f.  `oldEnergy` is updated with `testEnergy`.
  - **Arguments:**
    - `Optimizer :: EXTERNAL`: An external subroutine that optimizes the electron density.
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: The electron density in real space.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: Array of energy components; `energy(1)` is the total energy.
    - `forces :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output array for forces; `forces(:,:,1)` holds total forces.
    - `frozenIon :: LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN)`: Mask indicating frozen ions.

# Important Variables/Constants

- **Algorithm Control:**
    - `maxIonStep :: INTEGER`: (Imported from `IonOptimizers`) Maximum number of QuickMin iterations.
    - `timeStep :: REAL(KIND=DP)`: (Imported from `IonOptimizers`, but modified locally) The dynamic time step for the Verlet-like updates. It's adapted (increased or decreased) based on the outcome of each step.
    - `forceCutoff :: REAL(KIND=DP)`: (Imported from `IonOptimizers`) Convergence criterion; the optimization stops if the maximum force on any ion falls below this value.
- **State Variables:**
    - `velocity :: REAL(KIND=DP), DIMENSION(cell%numIon, 3)`: Current velocities of the ions in the QuickMin dynamics.
    - `lastCoord :: REAL(KIND=DP), DIMENSION(cell%numIon, 3)`: Stores the ion coordinates from the previous step, used to revert if an overstep occurs.
    - `minCoord :: REAL(KIND=DP), DIMENSION(cell%numIon, 3)`: Stores the coordinates corresponding to the state with the lowest maximum force encountered so far.
    - `initForces, prevForces :: REAL(KIND=DP), DIMENSION(SIZE(forces,1), SIZE(forces,2), SIZE(forces,3))`: Store forces from initial and previous steps.
    - `oldEnergy, testEnergy :: REAL(KIND=DP), DIMENSION(SIZE(energy))`: Store total energy from previous and current (trial) steps.
    - `overStep :: LOGICAL`: A flag indicating if the most recent step resulted in an energy increase or undesirable force-velocity alignment.
- **Parameters for Adaptive Timestep/Reset:**
    - The timestep is typically multiplied by 2 if a step is successful and divided by 5 if an overstep occurs.
    - Velocity is reset to zero if an overstep occurs.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision).
- **`TIMER`**: For performance timing (`TimerStart`, `TimerStop`).
- **`OUTPUT`, `OutputFiles`**: For writing messages and standard output (`WrtOut`, `outputUnit`).
- **`MPI_Functions`**: For `rankGlobal` (to control output from only rank 0).
- **`IonOptimizers`**: For configuration parameters (`maxIonStep`, initial `timeStep`, `forceCutoff`) and timer variables (`watch`, `watch2`).
- **`CellInfo`**: For `cell` derived type (to access `cell%numIon`, `cell%ionTable`, `cell%cellReal`).
- **`MathFunctions`**: For `Inverse` (to convert Cartesian steps to fractional for coordinate updates).
- **`Report`**: For formatted reporting of geometry optimization progress (`GeometryMinimizerReportHeader`, `GeometryMinimizerReportFooter`, `GeometryMinimizerReportSteps`).
- **`CalForces`**: For `CalculateForces` (to compute forces on ions at each step).
- **`RefreshIons`**: For `RefreshIonTerms` (called after ion positions are updated).

The `QuickMinOptimization` routine iteratively adjusts ion positions. In each step, it predicts new positions, calls the provided `Optimizer` to relax the electronic structure, calculates new forces, and then, based on the change in energy and the relationship between forces and velocities, decides how to adjust the `timeStep` and `velocity` for the next iteration. This makes it a self-adapting algorithm for finding local energy minima.
