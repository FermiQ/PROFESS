# Overview

The `NVE` module implements Molecular Dynamics (MD) simulations in the microcanonical ensemble. In the NVE ensemble, the Number of particles (N), the system Volume (V), and the total Energy (E) are conserved. This module uses the standard Velocity Verlet algorithm to integrate Newton's equations of motion.

The primary subroutine `RunNVE` orchestrates the simulation, including initialization of velocities, the main integration loop, and output of trajectory and thermodynamic data. It relies on utility functions from the `MolecularDynamics` module for many common tasks.

# Key Components

- **`MODULE NVE`**: The main container for NVE MD simulation logic.

- **`SUBROUTINE RunNVE(Optimizer, energy, forces)`**:
  - **Description:** This is the main driver routine for an NVE MD simulation.
    1.  Performs initialization: sets up timers, calculates total number of steps (`nStep`), allocates velocity arrays.
    2.  Initializes atomic velocities:
        - If `rstMD == 0` (new run): Calls `InitRandomSeed` and `InitVelocity` (from `MolecularDynamics`) to assign initial velocities based on the target `temperature`.
        - If `rstMD == 1` or `rstMD < 0` (restart): Calls `RestartMD` (from `MolecularDynamics`) to load velocities and positions.
    3.  Calculates initial forces by calling the external `Optimizer` (for electron density) and then `CalculateForces`.
    4.  Calculates and reports the initial conserved total energy. Initializes `cartNoWrap` for MSD.
    5.  Enters the main MD loop for `nStep` steps:
        a.  Optionally removes center of mass motion and adjusts center of mass position.
        b.  Calls `MonitorMeanSquareDisplacement`.
        c.  Performs the Velocity Verlet integration:
            i.  `velocity(t + dt/2) = velocity(t) + forces(t)/mass * dt/2`
            ii. `positions(t + dt) = positions(t) + velocity(t + dt/2) * dt`. Updates both unwrapped Cartesian (`cartNoWrap`) and wrapped fractional (`cell%ionTable%coord`) coordinates.
            iii. Calls `RefreshIonTerms`, then `Optimizer` (to get self-consistent density for new positions), then `CalculateForces` to get `forces(t + dt)`.
            iv. `velocity(t + dt) = velocity(t + dt/2) + forces(t + dt)/mass * dt/2`.
        d.  Calculates kinetic energy and the total conserved energy (`hamiltonian`).
        e.  Outputs trajectory data (`OutputMDGeom`, `OutputIonVelocityAndStates`).
        f.  Reports step information.
  - **Arguments:**
    - `Optimizer :: EXTERNAL`: An external subroutine to optimize the electron density at each MD step.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: Array of energy components from the electronic structure calculation. `energy(1)` is the total potential energy.
    - `forces :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output array for forces on ions.

- **`SUBROUTINE RescaleVelocity(KE, newKE, numIon)`** (Internal to `RunNVE`):
  - **Description:** Rescales all atomic velocities by a common factor `sqrt(newKE/KE)` to achieve a target total kinetic energy `newKE`.
  - **Arguments:**
    - `KE :: REAL(KIND=DP), INTENT(IN)`: Current total kinetic energy.
    - `newKE :: REAL(KIND=DP), INTENT(IN)`: Target total kinetic energy.
    - `numIon :: INTEGER, INTENT(IN)`: Number of ions.
  - **Note:** While present, the main NVE loop in `RunNVE` calculates the conserved energy but does not actively use this routine to enforce it by rescaling velocities each step; this would typically violate true NVE dynamics if `newKE` was fixed to maintain constant temperature.

- **`FUNCTION Conserved(KE, PE) RESULT(total_energy)`** (Internal to `RunNVE`):
  - **Description:** Calculates the total energy of the system, `E_total = KE + PE`, which should be conserved in an NVE simulation.
  - **Arguments:**
    - `KE :: REAL(KIND=DP), INTENT(IN)`: Total kinetic energy of the ions.
    - `PE :: REAL(KIND=DP), INTENT(IN)`: Total potential energy of the system (from `energy(1)`).
  - **Return Value:** `total_energy :: REAL(KIND=DP)`.

# Important Variables/Constants

This module primarily utilizes parameters and variables imported from the `MolecularDynamics` module, such as:
- `dt :: REAL(KIND=DP)`: The MD time step.
- `timeTot :: REAL(KIND=DP)`: The total simulation time.
- `temperature :: REAL(KIND=DP)`: Used for initializing velocities in a new run.
- `rstMD :: INTEGER`: Controls restart behavior.
- `startStep :: INTEGER`: Starting step number.
- `cartNoWrap`: Array for unwrapped Cartesian coordinates.
- `freedom`: Degrees of freedom.
- `watch, watch2`: Timer objects.

Key local variables in `RunNVE`:
- `velocity(:,:) :: REAL(KIND=DP), ALLOCATABLE`: Stores atomic velocities.
- `hamiltonian :: REAL(KIND=DP)`: Stores the calculated total conserved energy for the current step.
- `twiceKE :: REAL(KIND=DP)`: `2 * KE_atomic`, used in temperature calculation.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`MolecularDynamics`**: This is a foundational dependency, providing numerous parameters (like `dt`, `timeTot`, `temperature`, `rstMD`) and essential utility subroutines (e.g., `InitRandomSeed`, `InitVelocity`, `RestartMD`, `GetAtomKE`, `MonitorMeanSquareDisplacement`, `OutputMDGeom`, `OutputIonVelocityAndStates`, `RemoveMovementOfCenterOfMass`, `AdjustCenterOfMass`, `CalFreedom`).
- **`Constants`**: For `DP`, `fundamentalTime`, `boltzmann`.
- **`MathFunctions`**: For `Inverse` and `Vecmul` (used for coordinate transformations).
- **`Output`, `OutputFiles`**: For writing output (`WrtOut`, `outputUnit`).
- **`MPI_Functions`**: For `rankGlobal` (to control output) and `TITLE`.
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`, `watch`, `watch2`).
- **`CellInfo`**: For the `cell` derived type, providing access to cell and ion data.
- **`SYS`**: For `frozenIon` (mask for frozen atoms) and `rhoR` (electron density passed to `CalculateForces`).
- **`Report`**: For formatted reporting of MD progress (`GeometryMinimizerReportHeader`, `GeometryMinimizerReportSteps`, `GeometryMinimizerReportFooter`).
- **`CalForces`**: For `CalculateForces` to get forces on ions.
- **`RefreshCell`**: For `RefreshCellSetup` (called during `RestartMD` if cell parameters were part of the restart, though cell is fixed in NVE).
- **`RefreshIons`**: For `RescaleDensity` (called during `RestartMD`) and `RefreshIonTerms` (called after ion positions change).

The `RunNVE` subroutine is the main entry point. It requires an external `Optimizer` routine to update the electronic structure (and thus potential energy and forces) whenever ionic positions change. The NVE dynamics are propagated using the Velocity Verlet algorithm.
