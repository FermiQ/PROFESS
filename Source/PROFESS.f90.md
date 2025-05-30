# Overview

The `OFDFT` program is the main entry point and driver for the PROFESS (PRinceton Orbital-Free Electronic Structure Software) package. It controls the overall workflow of an orbital-free density functional theory calculation.

The program's responsibilities include:
1.  Initializing the computing environment (clocks, MPI for parallel runs).
2.  Reading the primary input (system name) from the command line.
3.  Calling routines to initialize all parameters, read detailed input files, and set up various computational modules.
4.  Printing a header with simulation details.
5.  Deciding whether to perform a static optimization (e.g., cell and/or ionic geometry relaxation) or a molecular dynamics run based on input parameters.
6.  Invoking the appropriate high-level routines for the chosen calculation type.
7.  Performing final cleanup and finalization of the MPI environment.

# Key Components

- **`PROGRAM OFDFT`**: The main program executable.

**Workflow:**
1.  **Initialization:**
    - Calls `InitClocks()` (presumably from `Timer` module) to set up timers.
    - Calls `StartClock('PROFESS')` to start timing the whole run.
    - Calls `InitializeMPI()` (from `MPI_Functions`) to set up the parallel environment and get `rankGlobal` and `sizeGlobal`.
    - Sets `outputRank` (from `Output` module, used for rank-specific printing) to `rankGlobal`.

2.  **Command-Line Argument Processing:**
    - MPI rank 0 reads the first command-line argument as `systemName`. This name is used as the base for input/output file names.
    - If `systemName` is not provided, the program stops.
    - `systemName` is broadcast to all other MPI processes using `BcastCharacter` (from `MPI_Functions`).

3.  **System Setup:**
    - Calls `InitPROFESS(systemName)` (from `Initializer` module). This comprehensive routine handles:
        - Reading all input option files (e.g., `systemName.inpt`).
        - Reading geometry files (e.g., `systemName.ion`).
        - Reading pseudopotential files.
        - Setting up various modules (FFT, KEDFs, cell parameters, etc.).
        - Reading or generating the initial electron density.
    - Calls `ReportHeader(systemName, sizeGlobal)` (from `REPORT` module) to print simulation setup information to the output file.

4.  **Main Calculation Stage:**
    - Checks the value of `mdType` (a variable likely set during `InitPROFESS` via input files, originating from the `Optimizer` module).
    - **If `mdType == -1` (default, indicating no molecular dynamics):**
        - Prints a message indicating static optimization.
        - Calls `OptimizeCell()` (from `Optimizer` module). This routine typically handles the relaxation of ionic positions and/or the simulation cell shape and volume to find the minimum energy structure.
    - **Else (if `mdType` indicates an MD run, e.g., NVE, NVT, NPT):**
        - Prints a message indicating a molecular dynamics run.
        - Calls `DoMolecularDynamics()` (from `Optimizer` module). This routine serves as a driver for the selected type of MD simulation.

5.  **Finalization:**
    - Calls `FinalizeOutput()` (from `REPORT` module) to print summary information, timings, etc.
    - Calls `CleanInitializeInputs()` (from `Initializer` module) to deallocate arrays and perform other cleanup from the initialization phase.
    - Calls `QuitMPI()` (from `MPI_Functions`) to shut down the MPI environment properly.
    - Calls `StopClock('PROFESS')` to record total execution time.

# Important Variables/Constants

- **`systemName :: CHARACTER(LEN=systemNameLen)`**: Stores the base name for the calculation, provided as a command-line argument. Used to construct various input and output file names.
- **`mdType :: INTEGER`**: (Imported from `Optimizer` module) An integer flag that determines the type of calculation:
    - `-1`: Static optimization (default).
    - Other values: Indicate different types of molecular dynamics ensembles (e.g., 0 for NVE, 1 for NVT, etc., specific values depend on `Optimizer` module's conventions).

# Usage Examples

To run the program, you would typically execute it from the command line, providing a system name:
```bash
./PROFESS.x my_system
```
This would expect input files like `my_system.inpt`, `my_system.ion`, etc.

# Dependencies and Interactions

The `PROGRAM OFDFT` links together many high-level modules:

- **`CONSTANTS`**: For `DP` (double precision) and `systemNameLen`.
- **`MPI_Functions`**: For MPI initialization (`InitializeMPI`), finalization (`QuitMPI`), rank/size information (`rankGlobal`, `sizeGlobal`), and broadcasting (`BcastCharacter`).
- **`OutputFiles`**: For `outputUnit` (used by `REPORT`).
- **`Output`**: For `outputRank` (though this is set using `rankGlobal` from `MPI_Functions`).
- **`REPORT`**: For printing header information (`ReportHeader`) and final summaries (`FinalizeOutput`).
- **`Initializer`**: For the main setup routine `InitPROFESS` and cleanup routine `CleanInitializeInputs`.
- **`Optimizer`**: For the high-level calculation drivers `OptimizeCell` (for static calculations) and `DoMolecularDynamics` (for MD), and for the `mdType` flag.
- **`Timer`** (Implicit): Calls to `InitClocks()`, `StartClock()`, and `StopClock()` suggest usage of a timing module, conventionally named `Timer` in this codebase.

The program follows a sequential flow: initialize, read inputs, perform main calculation (optimization or MD), and then clean up. The specific physics and algorithms are encapsulated within the modules called by `OFDFT`.
