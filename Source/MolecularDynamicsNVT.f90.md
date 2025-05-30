# Overview

The `NVT` module implements Molecular Dynamics (MD) simulations in the canonical ensemble, commonly referred to as NVT dynamics. In this ensemble, the Number of particles (N), the system Volume (V), and the average Temperature (T) are maintained constant.

To control the temperature, this module utilizes a Nose-Hoover thermostat, which is an extended Lagrangian method that introduces an additional degree of freedom (the thermostat variable `xLogS` and its conjugate momentum `vLogS` related to `Qmass`). The equations of motion for both the physical system and the thermostat variables are integrated using a reversible, explicit integrator, likely based on the work of Martyna et al., which often involves a multiple time-stepping approach (RESPA-like) via Yoshida-Suzuki decomposition for enhanced stability and accuracy.

# Key Components

- **`MODULE NVT`**: The main container for NVT MD simulation logic.

- **`SUBROUTINE RunNVT(Optimizer, energy, forces)`**:
  - **Description:** This is the primary driver subroutine for an NVT MD simulation.
    1.  **Initialization:** Sets up timers, calculates total number of MD steps (`nStep`), allocates arrays for `velocity` and `intCoeff` (integrator coefficients). Calls `MakeIntCoeff` to populate `intCoeff`.
    2.  **Velocity and Thermostat Initialization:**
        - If `rstMD == 0` (new run): Calls `InitRandomSeed` and `InitVelocity` to assign initial atomic velocities based on the target `temperature`. Thermostat variables `xLogS` and `vLogS` are initialized to zero.
        - If `rstMD == 1` or `rstMD < 0` (restart): Calls `RestartMD` (from `MolecularDynamics`, with `md_type=1`) to load atomic positions, cell, velocities, and thermostat variables (`xLogS`, `vLogS`) from previous run. Calls `RefreshCellSetup` and `RescaleDensity`.
    3.  Calculates initial forces (via `Optimizer` and `CalculateForces`) and the initial value of the conserved quantity for the Nose-Hoover Hamiltonian (`NHhamiltonian`).
    4.  Initializes `cartNoWrap` for Mean Square Displacement (MSD) calculations.
    5.  Prints initial energy and trajectory information.
    6.  **Main MD Loop:** Iterates for `nStep` steps:
        a.  Optionally removes center of mass motion (`RemoveMovementOfCenterOfMass`) and adjusts center of mass position (`AdjustCenterOfMass`).
        b.  Calls `MonitorMeanSquareDisplacement`.
        c.  **Integration Step (Velocity Verlet combined with Nose-Hoover):**
            i.  Calls `NHIntegrator` for the first half-update of thermostat variables and scaling of atomic velocities.
            ii. Updates atomic velocities for `dt/2`: `velocity(t + dt/2) = velocity_scaled(t) + forces(t)/mass * dt/2`.
            iii. Updates atomic positions for `dt`: `positions(t + dt) = positions(t) + velocity(t + dt/2) * dt`. Both unwrapped Cartesian (`cartNoWrap`) and wrapped fractional (`cell%ionTable%coord`) coordinates are updated.
            iv. Calls `RefreshIonTerms` (updates ion-ion energy, ion-electron potential).
            v.  Calls `Optimizer` (for electron density self-consistency).
            vi. Calls `CalculateStress` and `PrintStress`.
            vii. Calls `CalculateForces` to get `forces(t + dt)`.
            viii. Updates atomic velocities for the second `dt/2`: `velocity(t + dt) = velocity(t + dt/2) + forces(t + dt)/mass * dt/2`.
            ix. Calculates total atomic kinetic energy (`twiceKE = 2 * KE`).
            x.  Calls `NHIntegrator` again for the second half-update of thermostat variables and scaling of atomic velocities.
        d.  Outputs MD geometry (`OutputMDGeom`) and velocities/thermostat state (`OutputIonVelocityAndStates`).
        e.  Calculates and prints the conserved `NHhamiltonian` and other thermodynamic data.
        f.  Optionally reads a new target temperature via `ReadNewTemp`.
  - **Arguments:**
    - `Optimizer :: EXTERNAL`: External subroutine to optimize electron density.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: Energy components from electronic structure calculation.
    - `forces :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Forces on ions.

- **`SUBROUTINE NHIntegrator`** (Internal to `RunNVT`):
  - **Description:** Implements the integration steps for the Nose-Hoover thermostat variables (`xLogS`, `vLogS`) and applies the thermostatting effect to particle velocities. It uses a multiple time-step integration scheme (loops `iResn` times, and within each, `iYosh` times using `intCoeff`).
    In each sub-step:
    1.  Updates `vLogS` for a fraction of the time step using `gLogS = (twiceKE - freedom*temperature) / Qmass`.
    2.  Calculates a velocity scaling factor `scaling = scaling * EXP(- (sub_step_dt/2) * vLogS)`.
    3.  Updates `xLogS` using the current `vLogS`.
    4.  Recalculates `gLogS` with the (implicitly scaled) `twiceKE`.
    5.  Updates `vLogS` again for the second fraction of the sub-step.
    After all sub-steps, the final accumulated `scaling` factor is applied to all atomic velocities, and `twiceKE` is updated accordingly.
  - **Side Effects:** Modifies `xLogS`, `vLogS`, `velocity`, `twiceKE`.

- **`FUNCTION NHhamiltonian(KE, PE, xLogS, vLogS) RESULT(conserved_H)`** (Internal to `RunNVT`):
  - **Description:** Calculates the conserved quantity (Hamiltonian) for the Nose-Hoover NVT ensemble:
    `H_NH = KE_atomic + PE_system + (1/2)*Qmass*vLogS^2 + N_freedom*k_B*T*xLogS`.
  - **Arguments:**
    - `KE :: REAL(KIND=DP), INTENT(IN)`: Total kinetic energy of atoms.
    - `PE :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: Potential energy array; `PE(1)` is the total system potential energy.
    - `xLogS :: REAL(KIND=DP), INTENT(IN)`: Current "position" of the thermostat.
    - `vLogS :: REAL(KIND=DP), INTENT(IN)`: Current "velocity" of the thermostat.
  - **Return Value:** `conserved_H :: REAL(KIND=DP)`.

- **`LOGICAL FUNCTION check_isnan(a)`**:
  - **Description:** A utility function to check if a `REAL(KIND=DP)` value `a` is NaN (Not a Number). Returns `.TRUE.` if `a` is NaN, `.FALSE.` otherwise.
  - **Arguments:** `a :: REAL(KIND=DP), INTENT(IN)`.

# Important Variables/Constants

This module relies heavily on parameters and utility routines from the `MolecularDynamics` module, including:
- **MD Parameters:** `dt`, `timeTot`, `temperature`, `Qmass` (thermostat mass), `rstMD`, `nResn`, `nYosh`.
- **Atomic Data:** `velocity`, `cartNoWrap`.
- **System Info:** `freedom`.

Key local variables in `RunNVT`:
- `xLogS, vLogS :: REAL(KIND=DP)`: Thermostat position and velocity.
- `intCoeff(:) :: REAL(KIND=DP), ALLOCATABLE`: Coefficients for the Yoshida-Suzuki integrator.
- `twiceKE :: REAL(KIND=DP)`: Stores `2 * KE_atomic`.
- `hamiltonian :: REAL(KIND=DP)`: Stores the value of the conserved quantity.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`MolecularDynamics`**: This is a primary dependency, providing many basic MD parameters, data structures (like `cartNoWrap`), and utility subroutines (e.g., `InitVelocity`, `RestartMD`, `GetAtomKE`, `MakeIntCoeff`, `MonitorMeanSquareDisplacement`, `OutputMDGeom`, `OutputIonVelocityAndStates`, `RemoveMovementOfCenterOfMass`, `AdjustCenterOfMass`, `CalFreedom`).
- **`CellInfo`**: For the `cell` derived type.
- **`Constants`**: For physical constants (`DP`, `PI`, `fundamentalTime`, `boltzmann`, `hartreeToeV`).
- **`Output`, `OutputFiles`**: For writing output (`WrtOut`, `outputUnit`, `printStress`).
- **`MPI_Functions`**: For `rankGlobal` (to control output by a single processor) and `Title`.
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`, `stopwatch`).
- **`SYS`**: For `frozenIon` (mask for frozen atoms), `energy` array (from OFDFT), and `rhoR` (electron density).
- **`Report`**: For formatted reporting of MD progress (`GeometryMinimizerReportHeader`, `GeometryMinimizerReportSteps`, `GeometryMinimizerReportFooter`).
- **`RefreshCell`**: For `RefreshCellSetup` (used during `RestartMD`).
- **`CalStress`**: For `CalculateStress`.
- **`CalForces`**: For `CalculateForces`.
- **`RefreshIons`**: For `RefreshIonTerms` (called after particle moves) and `RescaleDensity` (used during `RestartMD`).

The `RunNVT` subroutine is the main entry point for NVT simulations. It requires an external `Optimizer` routine to update the electronic structure (and thus potential energy and forces) whenever ionic positions change during the MD steps. The NVT dynamics are propagated using a Velocity Verlet scheme for particles, tightly coupled with the Nose-Hoover thermostat integration performed by `NHIntegrator`.
