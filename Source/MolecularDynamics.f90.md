# Overview

The `MolecularDynamics` module provides a collection of utility subroutines that support the execution and analysis of Molecular Dynamics (MD) simulations. These routines handle tasks such as:
- Initializing simulations (velocities, random seeds).
- Restarting simulations from previous states.
- Managing and manipulating coordinates and velocities (e.g., handling center of mass motion, unwrapped coordinates for Mean Square Displacement).
- Calculating system properties like kinetic energy and degrees of freedom.
- Outputting trajectory data (positions, velocities, thermostat/barostat states).
- Dynamically adjusting simulation parameters like temperature.

This module serves as a foundational toolkit for more specific MD ensemble implementations (like NVE, NVT, NPT found in other modules).

# Key Components

- **`MODULE MolecularDynamics`**: The main container for MD utility routines.

- **`SUBROUTINE RestartMD(md_type, vel, xLogS, vLogS, vBoxG)`**:
  - **Description:** Reads data from restart files to continue a previous MD simulation. It reads ion positions and cell parameters from an "ion restart file" (e.g., `ion.restart` or `ion.step.dat`) and velocities, thermostat variables (`xLogS`, `vLogS`), and barostat variables (`vBoxG` for NPT) from a "velocity restart file" (e.g., `vel.restart` or `vel.step.dat`).
  - **Arguments:**
    - `md_type :: INTEGER, INTENT(IN)`: Type of MD ensemble (1:NVT, 2:NPT, 3:NVE) to guide reading.
    - `vel(3,cell%numIon) :: REAL(KIND=DP), INTENT(OUT)`: Atomic velocities.
    - `xLogS, vLogS :: REAL(KIND=DP), INTENT(OUT)`: Thermostat position and velocity.
    - `vBoxG(3,3) :: REAL(KIND=DP), OPTIONAL, INTENT(OUT)`: Barostat (cell velocity) variables for NPT.

- **`SUBROUTINE RemoveMovementOfCenterOfMass(vel, adjust)`**:
  - **Description:** Calculates the total momentum of the system and, if `adjust` is `.TRUE.`, subtracts the center of mass velocity from each atom's velocity to ensure zero net momentum.
  - **Arguments:**
    - `vel(3,cell%numIon) :: REAL(KIND=DP), INTENT(INOUT)`: Atomic velocities.
    - `adjust :: LOGICAL, INTENT(IN)`: If true, removes CM velocity.

- **`SUBROUTINE GetAtomKE(vel, ke)`**:
  - **Description:** Calculates the total classical kinetic energy of the atoms: `KE = SUM_i { 0.5 * m_i * |v_i|^2 }`.
  - **Arguments:**
    - `vel(3,size(cell%ionTable,1)) :: REAL(KIND=DP), INTENT(IN)`: Atomic velocities.
    - `ke :: REAL(KIND=DP), INTENT(OUT)`: Calculated total kinetic energy.

- **`SUBROUTINE MakeIntCoeff(nYosh, intCoeff)`**:
  - **Description:** Initializes coefficients for Yoshida-Suzuki reversible integrators of order `nYosh` (supports 1, 3, 5). These are used in some advanced MD integration schemes.
  - **Arguments:**
    - `nYosh :: INTEGER, INTENT(IN)`: Order of the integrator.
    - `intCoeff(nYosh) :: REAL(KIND=DP), INTENT(OUT)`: Array of integrator coefficients.

- **`SUBROUTINE InitRandomSeed()`**:
  - **Description:** Initializes the pseudo-random number generator seed using the system clock, ensuring different initial random numbers for different runs.

- **`SUBROUTINE InitVelocity(temperature, v)`**:
  - **Description:** Initializes atomic velocities for the start of an MD simulation. Velocities are drawn from a Maxwell-Boltzmann distribution corresponding to the target `temperature`. The routine then scales these velocities to precisely match the target temperature, after optionally removing center-of-mass motion and zeroing velocities of frozen atoms.
  - **Arguments:**
    - `temperature :: REAL(KIND=DP), INTENT(IN)`: Target temperature (in Hartree units).
    - `v(3,cell%numIon) :: REAL(KIND=DP), INTENT(OUT)`: Initialized atomic velocities.

- **`SUBROUTINE MonitorMeanSquareDisplacement(step, msd, diffuCoeff)`**:
  - **Description:** Calculates the Mean Square Displacement (MSD) of atoms relative to their positions at `msd_startTime`. Also computes the diffusion coefficient `D = MSD / (6*t)`. Uses unwrapped coordinates stored in `cartNoWrap`.
  - **Arguments:**
    - `step :: INTEGER, INTENT(IN)`: Current MD step.
    - `msd :: REAL(KIND=DP), INTENT(OUT)`: Calculated MSD (in Angstrom^2).
    - `diffuCoeff :: REAL(KIND=DP), INTENT(OUT)`: Calculated diffusion coefficient (in Angstrom^2/ps).

- **`SUBROUTINE OutputIonVelocityAndStates(iter, vel, md_type, xLogS, vLogS, vBoxG)`**:
  - **Description:** Writes the current atomic velocities, thermostat variables (`xLogS`, `vLogS`), and (for NPT) barostat variables (`vBoxG`) to a formatted data file (e.g., `vel.iter.dat`) for restarting or analysis. Output frequency is controlled by `dump_md_freq`.
  - **Arguments:** Similar to `RestartMD` but with `INTENT(IN)` for data.

- **`SUBROUTINE OutputMDGeom(iter)`**:
  - **Description:** Writes the current cell lattice vectors and unwrapped Cartesian atomic coordinates (`cartNoWrap`) to a formatted data file (e.g., `ion.iter.dat`). Also includes pseudopotential information. Output frequency is controlled by `dump_md_freq`.
  - **Arguments:** `iter :: INTEGER, INTENT(IN)`: Current MD iteration step.

- **`SUBROUTINE ReadNewTemp(temp, step)`**:
  - **Description:** If `fixTemperature > 0` and the current `step` is a multiple of `fixTemperature`, this routine reads a new target temperature from a file named `ChangeTemp.dat` and updates the simulation's target `temp`.
  - **Arguments:**
    - `temp :: REAL(KIND=DP), INTENT(INOUT)`: Target temperature, potentially updated.
    - `step :: INTEGER, INTENT(IN)`: Current MD step.

- **`SUBROUTINE AdjustCenterOfMass(shiftCM)`**:
  - **Description:** Calculates the center of mass of the system based on `cartNoWrap` coordinates. If `shiftCM` is `.TRUE.`, it shifts all atomic coordinates (both `cartNoWrap` and wrapped fractional `cell%ionTable%coord`) so that the center of mass is at the origin (0,0,0).
  - **Arguments:** `shiftCM :: LOGICAL, OPTIONAL, INTENT(IN)`.

- **`SUBROUTINE CalFreedom()`**:
  - **Description:** Calculates the total number of translational degrees of freedom (`freedom`) in the system, accounting for any constraints imposed by `frozenIon`.
  - **Output:** Updates the module-level variable `freedom`.

# Important Variables/Constants

- **MD Parameters (Module-Level):**
    - `tau_thermo, tau_baro :: REAL(KIND=DP)`: Characteristic fluctuation periods for thermostat and barostat couplings.
    - `timeTot :: REAL(KIND=DP)`: Total simulation time for the MD run.
    - `Qmass :: REAL(KIND=DP)`: "Mass" or inertia parameter for extended system variables (thermostats, barostats).
    - `dt :: REAL(KIND=DP)`: MD integration time step.
    - `temperature :: REAL(KIND=DP)`: Target temperature for NVT/NPT simulations (in Hartree units).
    - `nResn, nYosh :: INTEGER`: Parameters for Nose-Hoover chain integrators and Yoshida-Suzuki decomposition order.
- **Output & Restart Control (Module-Level):**
    - `output_stress_period :: INTEGER`: Frequency of stress tensor output.
    - `dump_md_freq :: INTEGER`: Frequency for dumping MD trajectory files (`ion.*.dat`, `vel.*.dat`).
    - `fixTemperature :: INTEGER`: Frequency for checking `ChangeTemp.dat` to update temperature.
    - `md_output_path :: CHARACTER(LEN=500)`: Directory path for MD output files.
    - `rstMD :: INTEGER`: Flag to control MD restart behavior.
    - `startStep :: INTEGER`: The initial step number of an MD run (can be > 0 if restarting).
- **MSD Related (Module-Level):**
    - `doMSDAtom :: LOGICAL`: If true, calculates and outputs per-atom MSD.
    - `msd_startTime :: REAL(KIND=DP)`: MD time (in steps) after which MSD calculation begins.
    - `msd, diffuCoeff :: REAL(KIND=DP)`: Calculated MSD and diffusion coefficient.
    - `cartNoWrap(:,:) :: REAL(KIND=DP), ALLOCATABLE`: Stores Cartesian coordinates of atoms without periodic wrapping, essential for correct MSD calculation.
    - `msd_coords0(:,:), msd_coordst(:,:) :: REAL(KIND=DP), ALLOCATABLE`: Store initial (at `msd_startTime`) and current unwrapped coordinates for MSD.
- **Other:**
    - `velRescale :: LOGICAL`: If true, initial velocities are rescaled to exactly match target temperature.
    - `freedom :: REAL(KIND=DP)`: Number of degrees of freedom.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For various physical constants (`fundamentalTime`, `PI`, `DP`, `BOHR`, `boltzmann`, `hartreeToeV`).
- **`MathFunctions`**: For `Inverse` (matrix inversion, e.g., for coordinate transformations) and `Vecmul`.
- **`Output`, `OutputFiles`**: For formatted output (`WrtOut`), error handling (`Error`), and file units (`outputUnit`, `errorUnit`).
- **`MPI_Functions`**: For `rankGlobal` (to control output by a single processor in parallel runs).
- **`CellInfo`**: For `Cell` derived type (accessing ion and cell data like `cell%ionTable`, `cell%elementTable`, `cell%cellReal`).
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`, `stopwatch` type).
- **`SYS`**: For `rhoR` (electron density, though not directly used by most utilities here) and `frozenIon` (mask for frozen atoms).

This module provides essential support functions for MD simulations. Specific MD algorithms (NVE, NVT, NPT) in other modules (`MolecularDynamicsNVE`, `MolecularDynamicsNPT`) will call these routines for initialization, I/O, and analysis tasks.
