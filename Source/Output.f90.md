# Overview

The `Output` module serves as a centralized hub for managing and generating all "pretty-printed" output files from the PROFESS (Princeton Orbital-Free Electronic Structure Software) application. It handles the formatting and writing of various types of simulation data, including energies, forces, stresses, system geometry (lattice and ion positions), electron densities, potentials, and Ewald summation details.

The module defines specific Fortran unit numbers for different output files and provides a set of logical flags to control what information is printed. It also includes utilities for MPI-aware output (ensuring only rank 0 writes certain messages or coordinating sequential writes for parallel data) and for gracefully terminating the program.

# Key Components

- **`MODULE Output`**: The main container for output routines and parameters.

- **Print Subroutines for Physical Quantities:**
    - `SUBROUTINE PrintForces(forces, systemName, ionStep)`: Prints atomic forces (in eV/A and Hartree/Bohr) to the main output unit or a separate `.force` file.
    - `SUBROUTINE PrintStress(stress)`: Prints the 3x3 stress tensor (in GPa and a.u.) and current lattice vectors to the main output.
    - `SUBROUTINE PrintEwald(...)`: Prints detailed information about Ewald summation parameters, energy components, and timings.
    - `SUBROUTINE PrintTransitionState()`: Outputs energy, forces, stress, and fractional ion coordinates to `transUnit` for interfacing with CINEB or similar transition state search codes.

- **Print Subroutines for System Geometry (typically to `outputGeomUnit`):**
    - `SUBROUTINE PrintLattice()`: Writes the `%BLOCK LATTICE_CART` section with cell vectors (in Angstroms).
    - `SUBROUTINE PrintIons(positions)`: Writes the `%BLOCK POSITIONS_FRAC` section with ion symbols and fractional coordinates. Can take optional `positions` argument.
    - `SUBROUTINE PrintPseudo()`: Writes the `%BLOCK SPECIES_POT` section listing element symbols and their pseudopotential filenames.
    - `SUBROUTINE PrintOptimizeIons()`: Writes the `%BLOCK ION_OPTIMIZATION` section with flags indicating which ionic degrees of freedom are optimized.
    - `SUBROUTINE PrintGeometry(cellRelaxFlag, Step, positions)`: A driver routine that calls `PrintLattice`, `PrintIons`, `PrintPseudo`, and `PrintOptimizeIons` to write a complete geometry file (e.g., `systemName.ion.step.geom` or `systemName.final.geom`).

- **Print Subroutines for Volumetric Data (Density, Potential):**
    - `SUBROUTINE Print3DArray(array, fileName, arrayUnit, recOffset, scaling)`: A generic, MPI-aware routine to write a 3D real array to a direct access formatted file. Handles sequential writing from multiple processors if data is distributed. Can scale down the grid for output if the file size exceeds `maxMB`.
    - `SUBROUTINE PrintDensity(density)`: Prints the 4D electron density `rhoR` (grid points x spin) to a file (e.g., `systemName.den`). Uses `Print3DArray` for each spin component.
    - `SUBROUTINE PrintPotential(potential)`: Prints a 4D potential (e.g., `potReal`) to a file (e.g., `systemName.pot`), similar to `PrintDensity`.
    - `SUBROUTINE PrintDensityClearly(rhoR, iter)`: Prints the electron density in a Tecplot-compatible formatted file. Handles parallel data gathering and writing.

- **Debugging and Utility Subroutines:**
    - **`INTERFACE Dump`**: Generic interface for dumping arrays for debugging purposes to `junkUnit`.
        - `SUBROUTINE DumpReal3D(array, systemName, increment)`
        - `SUBROUTINE DumpReal4D(array, systemName, increment)`
        - `SUBROUTINE DumpComplex(grid, systemName, preserve)` (dumps absolute values)
        - `SUBROUTINE DumpLogical(grid, systemName, preserve)`
    - `SUBROUTINE CleanOutput()`: Deallocates `outputIonTable`.
    - `SUBROUTINE WrtOut(unitOut, msg, nextLine)`: A general-purpose MPI-aware write routine. Only MPI rank 0 actually writes the message `msg` to the specified `unitOut`. `nextLine` controls advancing to a new line.
    - `SUBROUTINE QUIT(message)`: Terminates the program. Prints an optional `message`, finalizes MPI if active, and prints timing information via `printClock`.
    - `SUBROUTINE TITLE(message)` (Defined outside the module but in the same file): Prints a title string to standard output and the main output file from rank 0. (Currently disabled by an immediate `RETURN`).

# Important Variables/Constants

- **File Unit Parameters:**
    - `densityUnit, potentialUnit, outputGeomUnit, forceUnit, transUnit, junkUnit :: INTEGER, PARAMETER`.
- **Output Control Flags (Module-Level Booleans, default `.FALSE.`):**
    - `outputOptionGeom`: General flag for geometry output.
    - `outputFinalForces, outputEwald, outputFinalPotential, outputIntermediateDensity, outputFinalDensity, outputKernel, outputFinalGeometry, outputFinalStress, outputForcesSeparate, outputTransitionState`.
- **Output Formatting/Control:**
    - `maxMB :: REAL(KIND=DP)`: Maximum desired size in MB for density/potential files before down-scaling is applied (default: 2000 MB).
    - `geoOutputFreq, celOutputFreq :: INTEGER`: Frequency of geometry output during ionic or cell optimization steps respectively (default: -1, meaning off unless specifically set).
    - `outputSystemName :: CHARACTER(LEN=30)`: Base name used for constructing output filenames.
- **State Variables:**
    - `outputIonTable :: INTEGER, DIMENSION(:), POINTER`: Stores an unsorted list of ion indices to ensure consistent output order if the primary `ionTable` is sorted by type.
    - `outputRank :: INTEGER`: Stores the MPI rank of the current process (set from `MPI_Functions`).
- **Timing Arrays:**
    - `energyTime(:), singlePointTime(:), totalEnergyTime(:), totalPotentialTime(:) :: REAL(KIND=DP)`: Arrays to accumulate timing information for various parts of the calculation (likely populated by `Timer` module and used by `printClock` called in `QUIT`).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For physical constants (`hartreeToeV`, `bohr`, etc.) and `DP`.
- **`SYS`**: Provides access to core simulation data arrays that are to be printed, such as `rhoR` (real-space density), `potReal` (real-space potential), `forceIon`, `stress`, `energy`, and `frozenIon`.
- **`CellInfo`**: For `cell` derived type (providing lattice, ion positions, element types), and grid information like `m3G`, `n3G`, `n3Goff`.
- **`MATHFUNCTIONS`**: For `Volume`.
- **`TIMER`**: For `TimerStart`, `TimerStop`, `stopwatch` type, and `printClock` (used in `QUIT`).
- **`OUTPUTFILES`**: For global file unit numbers `outputUnit`, `errorUnit`, and the `outRank` variable.
- **`MPI_Functions`**: (Implicit or direct) For parallel-aware output. `Print3DArray` and `PrintDensityClearly` have explicit MPI logic. `WrtOut` and `QUIT` use `outputRank` (which is `rankGlobal` from `MPI_Functions`) to ensure rank 0 handles certain I/O.

This module is central to how users and developers observe the state and results of a PROFESS simulation. Its routines are called from various parts of the code, especially after significant computational steps (e.g., SCF convergence, MD step, geometry optimization step) or at the end of a run.
