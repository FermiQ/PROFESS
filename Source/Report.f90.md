# Overview

The `Report` module is dedicated to generating formatted, human-readable output that summarizes the parameters, progress, and results of an OFDFT (PROFESS) simulation. It provides a suite of subroutines to print various sections of the output file, including:
- A main header with simulation setup details.
- Step-by-step progress for electronic density minimization loops.
- Step-by-step progress for ionic geometry optimization or molecular dynamics loops.
- Detailed breakdowns of energy components and their calculation times.
- Final energies, forces, stresses, and geometry.
- Overall timing information for different parts of the calculation.

The output routines are generally MPI-aware, with MPI rank 0 typically responsible for writing most of the information to ensure clean and non-redundant output in parallel runs.

# Key Components

- **`MODULE Report`**: The main container for all reporting subroutines.

- **`SUBROUTINE ReportHeader(systemName, numProc)`**:
  - **Description:** Prints the main header at the beginning of the simulation output. This includes the program title, start date/time, system name, cell information (number of ions, types, lattice vectors, boundary conditions), and core simulation parameters (kinetic energy cutoff, number of spins, total electrons, grid dimensions, number of processors, and chosen algorithms for electronic and ionic minimization/dynamics).
  - **Arguments:**
    - `systemName :: CHARACTER(LEN=*), INTENT(IN)`: The base name for the simulation.
    - `numProc :: INTEGER, INTENT(IN)`: Total number of MPI processors.

- **`SUBROUTINE ReportAnswers(energy, title)`**:
  - **Description:** Prints a detailed breakdown of energy components (Thomas-Fermi, von Weizsacker, Wang-Teter/WGC/LQ/HQ, Exchange-Correlation, Coulombic, Ion-Electron, Ion-Ion) along with the time taken to calculate each component. It also prints total kinetic, total potential, and total energy.
  - **Arguments:**
    - `energy :: REAL(kind=DP), DIMENSION(:), INTENT(IN)`: Array of energy components (from `SYS` module).
    - `title :: CHARACTER(len=*), INTENT(IN)`: A title string for this energy report section.

- **`SUBROUTINE MinimizerReportHeader`**:
  - **Description:** Prints the header section for an electronic density minimization loop, indicating the type of algorithm being used (e.g., Sqrt Conjugate Gradient, Sqrt Truncated Newton).
  - **Side Effects:** Initializes `timeCumul` for the minimizer loop.

- **`SUBROUTINE MinimizerReportSteps(step, energy, potentialNorm, numEnergyLine, numEnergyBrack, restart, success, duration)`**:
  - **Description:** Prints a line of information for each step of an electronic density minimization loop. Includes step number, total energy, norm of the potential (gradient), number of LCG steps/line search evaluations, line search step size, cumulative FFT count, and cumulative time. Notes restarts (`+`) or problematic line searches (`*`, `^`). Optionally calls `ReportAnswers` and `PrintDensity` if verbosity levels are high.
  - **Arguments:** Various parameters detailing the progress of the current minimization step.

- **`SUBROUTINE MinimizerReportFooter(steps, cause, energy, numEnergyEvals, numPotentialEvals, duration)`**:
  - **Description:** Prints a summary at the conclusion of an electronic density minimization loop. Includes the reason for stopping (e.g., convergence, max steps exceeded), total steps, total energy/potential evaluations, machine precision for energy, and total minimization time. May call `ReportAnswers` if verbosity is high or if the minimizer failed.

- **`SUBROUTINE GeometryMinimizerReportHeader(extraInfo)`**:
  - **Description:** Prints the header section for an ionic geometry optimization or molecular dynamics loop, indicating the method used (Quickmin, CG, BFGS, NVT, NPT, NVE). Can include optional `extraInfo`.

- **`SUBROUTINE GeometryMinimizerReportSteps(step, duration, forces, maxForce, stepSize, overStep, velocityMag, stepDone, positions)`**:
  - **Description:** Prints a line of information for each step of an ion geometry optimization or MD run. Content varies slightly based on whether it's optimization or MD. For optimization, includes step number, max force, step size, velocity magnitude (if applicable), overstep flag, and time. For MD, includes step number, max force, MD time step, Hamiltonian/conserved quantity, temperature, and time. Optionally calls `PrintForces` and `PrintGeometry` based on verbosity and output frequency settings.
  - **Arguments:** Various parameters detailing the progress of the current geometry/MD step.

- **`SUBROUTINE GeometryMinimizerReportFooter(duration)`**:
  - **Description:** Prints a summary at the conclusion of an ion geometry optimization or MD loop, primarily the total time taken.

- **`SUBROUTINE FinalizeOutput()`**:
  - **Description:** Performs the final output tasks at the end of the entire PROFESS run.
    1.  Prints the final electron density (calls `PrintDensity` from `Output` module), adding back core density if density decomposition was used.
    2.  Calls `ReportAnswers` to print final energies.
    3.  Prints final forces (calls `PrintForces` from `Output`) if requested and not already printed.
    4.  Prints final stress (calls `PrintStress` from `Output`) if requested.
    5.  Prints final geometry (calls `PrintGeometry` from `Output`) if requested.
    6.  Stops the main 'PROFESS' timer.
    7.  Calls `PrintClock` (from `Timer` module) to print a summary of all timers.
    8.  Prints a run completion message with date and time.
    9.  Closes `outputUnit`, `errorUnit`, and `outputGeomUnit`.
    10. If `outputTransitionState` is true, calls `PrintTransitionState` (from `Output`).
    11. Performs an MPI barrier.

# Important Variables/Constants

This module primarily formats and prints data. Most of the important data variables are imported from other modules, particularly:
- `energy(:)` from `SYS`: For energy breakdowns.
- `forces(:,:,:)`, `stress(:,:)` from `SYS`: For final forces and stress.
- `cell` from `CellInfo`: For lattice and ion information.
- `outputRank`, `outputSystemName`, `lineLength`, various output control flags (e.g., `outputMinimizeDensity`, `geoOutputFreq`), and file units from `Output` and `OutputFiles`.
- `energyTime(:)` from `Output`: For timing of energy components.
- `kinetic` from `CellInfo`, `exchangeCorrelationFunct` from `Output` (originally `XC_LDA`): To label energy components correctly.
- `outputRhoMethod`, `outputIonMethod` from `Optimizer`: To label algorithm types in headers.
- `iCountFFT` from `Fourier`: For FFT statistics.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`OutputFiles`**: For `outputUnit`.
- **`MPI_Functions`**: For `rankGlobal`.
- **`OUTPUT`**: This module is heavily coupled with `OUTPUT`, relying on it for many boolean flags that control printing verbosity and frequency (e.g., `outputMinimizeDensity`, `outputFinalForces`), file unit definitions, `outputRank`, `outputSystemName`, `lineLength`, and even for some printing subroutines like `PrintDensity`, `PrintForces`, `PrintStress`, `PrintGeometry`, `PrintTransitionState`. It also uses timing arrays like `energyTime` from `OUTPUT`.
- **`CellInfo`**: For `cell` derived type (to access lattice, ion data) and `m3G` (global Z grid dimension).
- **`CONSTANTS`**: For `DP` and various conversion factors (`hartreeToeV`, `bohr`, `boltzmann`, `auToGPa`).
- **`SYS`**: For accessing the final state variables like `energy`, `rhoR`, `forceIon`, `stress` during `FinalizeOutput`.
- **`PlaneWave`**: For `energyCutoff` (used in `ReportHeader`).
- **`KEDF_DenDec`**: For `do_den_dec` flag and `AddCoreDensity` subroutine, used in `FinalizeOutput` when printing the final density.
- **`TIMER`**: For `TimerStart`, `TimerStop`, `stopwatch` type, and `PrintClock`, `PrintClockWith` (used in `FinalizeOutput`).
- **`Optimizer`**: For `outputRhoMethod` and `outputIonMethod` flags to describe the algorithm in headers.
- **`Fourier`**: For `iCountFFT` to report cumulative FFT counts.

The `Report` module's subroutines are called at different stages of a PROFESS run:
- `ReportHeader`: Called once near the beginning by the main program.
- `MinimizerReportHeader/Steps/Footer`: Called by electronic structure optimization routines.
- `GeometryMinimizerReportHeader/Steps/Footer`: Called by ionic geometry optimization or molecular dynamics routines.
- `ReportAnswers`: Called by minimizer footers or `FinalizeOutput` to detail energy components.
- `FinalizeOutput`: Called once at the very end of the main program.
