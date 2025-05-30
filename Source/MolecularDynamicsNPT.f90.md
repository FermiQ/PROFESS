# Overview

The `NPT` module implements Molecular Dynamics (MD) simulations in the isothermal-isobaric (NPT) ensemble, where the Number of particles (N), system Pressure (P), and system Temperature (T) are kept constant on average. This is achieved using an extended system approach, coupling the physical system to a Nose-Hoover thermostat (to control temperature) and a barostat (Parrinello-Rahman like, to control pressure by allowing the simulation cell volume and shape to fluctuate).

The integration of the equations of motion is based on the explicit reversible integrators for extended system dynamics developed by Martyna et al. (Mol. Phys. 87, 1117 (1996)). This often involves a multiple time-step algorithm (RESPA-like) where different parts of the system (e.g., atomic motion, thermostat evolution, barostat evolution) are updated on different time scales, typically using Yoshida-Suzuki decomposition for the propagators.

# Key Components

- **`MODULE NPT`**: The main container for NPT MD simulation logic.

- **`SUBROUTINE NPTFullCell(RhoOptimizer)`**:
  - **Description:** This is the main driver subroutine for an NPT MD simulation.
    1.  Initializes parameters, degrees of freedom, integrator coefficients, and random number seed.
    2.  Initializes atomic velocities and extended system variables (thermostat position `xLogS`, velocity `vLogS`; barostat/cell velocities `vBoxG`), either from random distributions (for a new run) or by reading from restart files (via `RestartMD` from `MolecularDynamics` module).
    3.  Ensures the center of mass is at the origin and applies box constraints.
    4.  Performs an initial calculation of energy, forces, stress, and the conserved quantity for the NPT ensemble.
    5.  Enters the main MD loop for `nStep` steps:
        a.  Calls `NPTFullCell_Integrator` for the first half-step evolution of velocities and extended variables.
        b.  Calls `MVV` (Modified Velocity Verlet) to update atomic positions and cell parameters, and perform a full update of atomic velocities.
        c.  Calls `NPTFullCell_Integrator` again for the second half-step evolution.
        d.  Calculates kinetic energies, the conserved quantity, and outputs trajectory/thermodynamic data.
  - **Arguments:**
    - `RhoOptimizer :: EXTERNAL`: An external subroutine to optimize electron density at each MD step.

- **`SUBROUTINE MVV(RhoOptimizer, dt, vel, vBoxG)`**:
  - **Description:** Implements a Modified Velocity Verlet-like step for updating particle positions, cell shape, and particle velocities.
    1.  Updates atomic velocities for the first half step: `vel(0) -> vel(dt/2)` using current forces.
    2.  Updates atomic Cartesian positions (`cartNoWrap`) and the cell matrix (`cell%cellReal`) based on current atomic velocities and cell velocities `vBoxG`, using an exponential propagator derived from diagonalizing `vBoxG`.
    3.  Recalculates fractional coordinates, updates cell-dependent quantities (`RefreshCellSetup`), and rescales density (`RescaleDensity`).
    4.  Calls `RhoOptimizer` to get the new electron density.
    5.  Calculates new atomic forces using `CalculateForces`.
    6.  Updates atomic velocities for the second half step: `vel(dt/2) -> vel(dt)` using the new forces.
  - **Arguments:** `RhoOptimizer`, `dt`, `vel` (INOUT), `vBoxG` (IN).

- **`SUBROUTINE NPTFullCell_Integrator(energy, intCoeff, xLogS, vLogS, vBoxG, velocity, qtStress)`**:
  - **Description:** Integrates the equations of motion for the extended system variables (thermostat, barostat) and atomic velocities for a half time step, using a multiple time-stepping scheme based on `nResn` (number of RESPA steps) and `nYosh` (order of Yoshida-Suzuki decomposition) with coefficients `intCoeff`. It updates `vLogS`, `xLogS`, `vBoxG`, and `velocity`.
  - **Arguments:** `energy` (IN), `intCoeff` (IN), `xLogS` (INOUT), `vLogS` (INOUT), `vBoxG` (INOUT), `velocity` (INOUT), `qtStress` (OUT, OFDFT stress).

- **`SUBROUTINE FrozeAtom(dataIn, flag)`**:
  - **Description:** Sets components of `dataIn` to zero for atoms that are marked as frozen in the `frozenIon` array. `flag=1` for velocities, `flag=2` for forces.
  - **Arguments:** `dataIn(:,:)` (INOUT), `flag` (IN).

- **`SUBROUTINE NPTFullCell_Conserved(PE, xLogS, vLogS, d)`**:
  - **Description:** Calculates the conserved quantity for the NPT ensemble (Martyna et al. formalism). This includes atomic kinetic and potential energy, `P_ext * Volume` term, barostat kinetic energy, and thermostat kinetic and potential energies. Also calculates and prints the effective temperature.
  - **Arguments:** `PE` (IN, potential energies), `xLogS` (IN), `vLogS` (IN), `d` (IN, system dimensionality = 3).

- **`SUBROUTINE GetAtomKETensor(vel, keTensor)`**:
  - **Description:** Calculates the atomic kinetic energy tensor: `KE_ab = SUM_i { 0.5 * m_i * v_i_a * v_i_b }`.
  - **Arguments:** `vel` (IN), `keTensor(3,3)` (OUT).

- **`SUBROUTINE GetBoxKE(bMass, vBoxG, boxKE)`**:
  - **Description:** Calculates the kinetic energy of the barostat (cell degrees of freedom): `KE_box = 0.5 * bMass * Trace(vBoxG^T * vBoxG)`.
  - **Arguments:** `bMass` (IN), `vBoxG` (IN), `boxKE` (OUT).

- **`SUBROUTINE BoxConstraint(force_box, constr_type)`**:
  - **Description:** Applies constraints to the forces (or velocities) acting on the cell's degrees of freedom based on `constr_type`. For example, `constr_type=1` keeps the box orthogonal by zeroing off-diagonal elements of `force_box`. Other types can fix specific dimensions.
  - **Arguments:** `force_box(3,3)` (INOUT), `constr_type` (IN).

- **`SUBROUTINE EstimateQMass(Nf, temp, tau, Qmass)`**:
  - **Description:** Estimates the thermostat mass `Qmass = Nf * kT * tau_thermo^2`.
  - **Arguments:** `Nf` (degrees of freedom), `temp` (temperature in Hartree), `tau` (thermostat period in a.u.), `Qmass` (OUT).

- **`SUBROUTINE EstimateBMass(Nf, temp, tau, Bmass)`**:
  - **Description:** Estimates the barostat mass `Bmass = (Nf + d) * kT * tau_baro^2`.
  - **Arguments:** `Nf`, `temp`, `tau` (barostat period in a.u.), `Bmass` (OUT).

- **`SUBROUTINE CheckNPTFullCellParam()`**:
  - **Description:** Checks if all necessary NPT parameters (e.g., `extPres`, `dt`, `temperature`, masses) are initialized. If `tau_thermo` or `tau_baro` are set, it calls `EstimateQMass` or `EstimateBMass` to set the corresponding masses. Prints a summary of NPT parameters.

- **`SUBROUTINE OutputStress(iter, stress_ofdft, stress_ke)`**:
  - **Description:** Outputs the OFDFT stress tensor, the kinetic stress tensor, and the total stress tensor (sum of both) to the main output file. Frequency controlled by `output_stress_period`.
  - **Arguments:** `iter` (IN), `stress_ofdft(3,3)` (IN), `stress_ke(3,3)` (IN).

- **`SUBROUTINE StressTensorKE(vel, keStress)`**:
  - **Description:** Calculates the kinetic contribution to the stress tensor: `sigma_ke_ij = (1/Volume) * SUM_k { m_k * v_k(i) * v_k(j) }`.
  - **Arguments:** `vel` (IN), `keStress(3,3)` (OUT).

# Important Variables/Constants

- **Module-Level NPT Parameters:**
    - `extPres :: REAL(KIND=DP)`: Target external pressure (in Hartree/Bohr^3).
    - `bMass :: REAL(KIND=DP)`: Mass parameter for the barostat. Can be estimated from `tau_baro`.
    - `constr_type :: INTEGER`: Flag specifying constraints on cell deformation (e.g., -1 for none, 1 for orthogonal).
- **Key Variables from `MolecularDynamics` Module:**
    - `tau_thermo, tau_baro :: REAL(KIND=DP)`: Characteristic time constants for thermostat and barostat.
    - `timeTot :: REAL(KIND=DP)`: Total simulation time.
    - `Qmass :: REAL(KIND=DP)`: Mass parameter for the thermostat. Can be estimated from `tau_thermo`.
    - `dt :: REAL(KIND=DP)`: MD time step.
    - `temperature :: REAL(KIND=DP)`: Target temperature (in Hartree).
    - `nResn, nYosh :: INTEGER`: Parameters for the multiple time-step integrator.
    - `freedom :: REAL(KIND=DP)`: Number of degrees of freedom.
    - `cartNoWrap`: Unwrapped Cartesian coordinates for MSD.
- **Extended System Variables (Local to `NPTFullCell`):**
    - `xLogS :: REAL(KIND=DP)`: Position of the Nose-Hoover thermostat.
    - `vLogS :: REAL(KIND=DP)`: Velocity of the Nose-Hoover thermostat.
    - `vBoxG(3,3) :: REAL(KIND=DP)`: Velocities of the cell lattice vectors (barostat variables).
- **Thermodynamic Quantities:**
    - `effTemp :: REAL(KIND=DP)`: Calculated effective temperature of the system.
    - `atomKE :: REAL(KIND=DP)`: Total kinetic energy of atoms.
    - `boxKE :: REAL(KIND=DP)`: Kinetic energy of the barostat.
    - `conserved :: REAL(KIND=DP)`: The conserved quantity in the extended NPT Hamiltonian.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`MolecularDynamics`**: This is a primary dependency, providing many basic MD parameters (like `dt`, `temperature`, `Qmass`, `rstMD`) and utility subroutines (e.g., `RestartMD`, `GetAtomKE`, `MakeIntCoeff`, `InitVelocity`, `OutputIonVelocityAndStates`, `OutputMDGeom`, `AdjustCenterOfMass`, `CalFreedom`).
- **`CellInfo`**: For the `cell` derived type, providing access to cell and ion information.
- **`Constants`**: For physical constants (`fundamentalTime`, `PI`, `DP`, `BOHR`, `boltzmann`, `hartreeToeV`).
- **`MathFunctions`**: For matrix operations like `Inverse`, `Vecmul`, `Symmetrize`, `Diag`, and `Det`.
- **`Output`, `OutputFiles`**: For writing output (`WrtOut`, `outputUnit`, `message`, `Error`).
- **`MPI_Functions`**: For `rankGlobal` to manage rank-specific operations like file I/O.
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`, `stopwatch`).
- **`SYS`**: For global system variables like `energy` (total energy components from OFDFT), `rhoR` (electron density), `forceIon` (forces from OFDFT), `frozenIon`.
- **`CalStress`**: For `CalculateStress` to get the OFDFT contribution to the virial.
- **`RefreshCell`**: For `RefreshCellSetup` to update cell-dependent quantities after cell changes.
- **`RefreshIons`**: For `RescaleDensity` to ensure correct electron count after cell volume changes.

The `NPTFullCell` routine orchestrates the simulation, calling specialized subroutines for different parts of the Martyna-Tuckerman-Tobias-Klein (MTTK) or similar NPT integration scheme. It relies on an external `RhoOptimizer` to update the electronic structure at each MD step.
