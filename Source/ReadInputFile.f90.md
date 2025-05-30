# Overview

The `ReadInputFile` module is responsible for parsing the main control input file for a PROFESS simulation, typically named with a `.inpt` extension (e.g., `systemname.inpt`). This file contains various keywords that specify the parameters and options for the calculation.

The module defines two main subroutines:
- `ReadOptions`: Reads through the input file, identifies keywords, and assigns the corresponding values to global variables located in many other modules throughout the PROFESS codebase. This is the primary mechanism for user control over the simulation.
- `CheckOptions`: Performs basic validation of some of the critical parameters read by `ReadOptions` to ensure they are within sensible ranges or that conflicting options have not been set.

This module acts as a central hub for configuring the entire simulation by setting parameters across a wide array of other specialized modules.

# Key Components

- **`MODULE ReadInputFile`**: The main container module.

- **Module-Level Parameters & Variables (for file naming):**
    - `defaultPseudo, defaultTrans, defaultInput, defaultGeometry, defaultDensity, defaultOutput, defaultError :: CHARACTER, PARAMETER`: Define default file extensions.
    - `transFile, inputFile, geometryFile, densityFile, outputFile, errorFile :: CHARACTER(...)`: Store the full names of various input/output files, constructed from the base `systemName` (provided externally, typically from command line) and the default extensions.
    - `fileDensityFlag :: LOGICAL`: Set to `.TRUE.` if a density file is specified via `RHOFILE` or `OLDDENS`.
    - `fileKernelFlag :: LOGICAL`: (Declared but not set by any keyword in the provided code).
    - `atomicDensityFlag :: LOGICAL`: Set to `.TRUE.` if `RHOATOM` keyword is present.
    - `xRepeat, yRepeat, zRepeat :: INTEGER`: Repetition factors for reading density from a file (default: 1).
    - `output_all_files :: INTEGER`: Flag for parallel output (0 for rank 0 only, 1 for all ranks; default 0).

- **`SUBROUTINE ReadOptions`**:
  - **Description:** This is the core routine that parses the `.inpt` file. It reads the file line by line. For each line, it attempts to match the first word (keyword) and then reads the subsequent option or value(s). It uses a large `SELECT CASE` structure to handle numerous keywords.
  - **Keywords Handled (examples):**
    - **General Setup:** `ECUT` (energy cutoff), `GDEN` (grid density), `DIME` (FFT dimension padding type).
    - **Task Control:** `MINIMIZE` (with options `RHO`, `IONS`, `CELL`), `CALCULATE` (with `FORCES`, `STRESS`).
    - **Algorithm Selection:** `METHOD` (for `RHO`, `IONS`, `CELL`), `DENCGALG` (CG type for density), `IONCGALG` (CG type for ions), `MBFGS`.
    - **Convergence Criteria:** `TOLENERGY`, `TOLPOTENTIAL`, `TOLSTRESS`, `TOLFORCE`, `TOLMAGMOMENT`.
    - **Iteration Limits:** `MAXITERDEN`, `MAXITERION`, `MAXITERCEL`, `MAXITERELEC`, `MAXITERMAG`.
    - **MD Parameters:** `MDOPATH`, `RSTMD`, `NVEMD` (TEMP, TIMETOT, DTIMESTEP), `NOSEHOOVER` (TEMP, TIMETOT, DTIMESTEP, QMASS, NRESN, NYOSH), `NPTMD` (TEMP, TIMETOT, DTIMESTEP, QMASS, BMASS, EXTPRESS, NRESN, NYOSH, TAUT, TAUB, CONSTRTYPE), `VELRESCALE`, `MSDINIT`, `MSDATOM`.
    - **KEDF Selection & Parameters:** `KEDFINETIC`, `PARAMKEDF` (LAMBDA, MU, ALPHA, BETA, GAMMA, CATA, CATB, CATG, RHO0, RHOS, MRHOS, LUMGPEXP, LUMGPFAC, PBECUTOFF, ALPHAWGC5, BETAWGC5, HCLAMBDA, CUTOFFRHO, RATIOCAT, ATF, BVW, MSHIFT, SCFC, RHOCUT, KLOC, ALOC, TOLK, CPENALTY, VWMODEL, LKT, LMUMGP, LBETMGP, LLAMMGP, LSIGMGP).
    - **XC Functional:** `EXCHCORR` (LDA, PBE).
    - **Ewald Sum Parameters:** `EWALDTOL`, `EWALDMAXREAL`, `EWALDETAINC`, `EWALDRECIPINC`.
    - **Spline Options:** `SPLINE` (ALL, IE, II, NONE, with order).
    - **Density Input:** `RHOFILE`, `OLDDENS`, `RHOATOM`, `RHOREPEAT` (X, Y, Z).
    - **Output Control:** `PRINT` (MINRHO, MINION, CELLGEOMFREQ, EWALD, STRESS, RESTARTINFO, RHOFINAL, IONFINAL, FORCES, MDINFO, KERNEL), `TRANSITION`.
    - **Misc./Debugging:** `RHOMAXMB`, `HOLDVAL` (RHO0, RHOS), `RHOVAC`, `RHOSTEPVAC`, `CHEATOPT`, `CHGTEMP`, `ALLFILES`.
  - **Side Effects:** Modifies a vast number of global parameters stored in various `USE`d modules.

- **`SUBROUTINE CheckOptions`**:
  - **Description:** Performs basic validation on some of the parameters read by `ReadOptions`. For example, it checks if `energyCutoff` or `gridSpacing` is set, if MD parameters are positive, and if `rho0` and `rhoS` are set when `bvac` is true. If an invalid or missing critical parameter is found, it calls the `Error` routine (from `OutputFiles`) to print messages and terminate the program.

# Important Variables/Constants

This module acts as the primary parser for user inputs. The importance of variables lies in how they control the behavior of the entire PROFESS suite. See the "Keywords Handled" section above for a list of what kind of parameters are set. The actual variables are typically located in the modules they configure (e.g., `energyCutoff` in `PlaneWave`, `kinetic` in `CellInfo`, `rhoMethod` in `Optimizer`, etc.).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

This module has **extensive dependencies** as it is responsible for setting parameters across almost all other functional modules in the PROFESS code.
- **`CONSTANTS`**: For `DP`, `PI`, unit conversion factors (`auPressure`, `boltzmann`, etc.), `systemNameLen`.
- **`OutputFiles`**: For `outputUnit`, `errorUnit`, `Error` interface, `message` buffer.
- **`Output`**: For `WrtOut` and numerous output control flags it sets (e.g., `outputFinalDensity`, `outputMinimizeGeometry`).
- **`ReadIonFile`**: For `inputUnit` parameter.
- **All configuration-receiving modules**: This includes `CellInfo`, `Sys`, `PlaneWave`, `SetupFFT`, `Ewald`, `IonElectronSpline`, `NearestDistance`, `RefreshIons`, `CBSpline`, all `KEDF_*` modules, `Optimizer`, `RhoOptimizers`, `RhoDirCG`, `RhoDirBFGS`, `IonOptimizers`, `CellOptimizers`, `MolecularDynamics`, `NPT`. The `USE` statements at the beginning of the module list these explicitly.
- **`MathFunctions`**: For `UpperCase`.
- **`MPI_Functions`**: For `rankGlobal` (to control which process reads/writes certain messages).

The `ReadOptions` subroutine is called early in the program initialization (by `Initializer.InitPROFESS`) immediately after the base `systemName` is determined. `CheckOptions` is called after `ReadOptions` to validate the settings. The parameters read by this module dictate the entire subsequent workflow and algorithmic choices of the simulation.
