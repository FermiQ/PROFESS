# Overview

The `Optimizer` module serves as the primary high-level controller for the type of calculation performed within the PROFESS (Princeton Orbital-Free Electronic Structure Software) package. It acts as a dispatcher, deciding whether to perform static calculations (like electron density minimization, ionic position optimization, or cell relaxation) or to run a molecular dynamics (MD) simulation.

The choice of action is determined by module-level parameters, primarily `mdType` and the `calOption` array, which are typically set based on user input read during the initialization phase. This module then calls the appropriate specialized routines from other modules that implement the specific algorithms.

# Key Components

- **`MODULE Optimizer`**: The main container module.

- **Module-Level Control Parameters:**
    - **`calOption(5) :: LOGICAL`**: An array of boolean flags that specify which tasks to perform.
        - `calOption(1)`: Minimize the electron density (default: `.TRUE.`).
        - `calOption(2)`: Minimize ionic positions inside the unit cell (default: `.FALSE.`).
        - `calOption(3)`: Relax the cell by minimizing the stress tensor (default: `.FALSE.`).
        - `calOption(4)`: Calculate the forces on the ions (default: `.FALSE.`).
        - `calOption(5)`: Calculate the stress tensor (default: `.FALSE.`).
    - **`cellRelax :: INTEGER`**: Specifies dimensions for cell relaxation if `calOption(3)` is true (default: -1, meaning full cell relaxation or as per `SteepestDecent`'s default).
        - `1`: Optimize only x-direction.
        - `2`: Optimize only y-direction.
        - `3`: Optimize only z-direction.
    - **`ionMethod :: INTEGER`**: Selects the algorithm for ion position optimization if `calOption(2)` is true (default: -1).
        - `0`: No minimization.
        - `2`: Quickmin (`IonOptQui.QuickMinOptimization`).
        - `3`: Legacy Conjugate Gradient (`IonOptCG.GradientOptimization`).
        - `4`: Conjugate Gradient version 2 (`IonOptCG2.ConjugateGradient`).
        - `5`: BFGS (`IonOptBFGS.IonLBFGS`).
    - **`rhoMethod :: INTEGER`**: Selects the algorithm for electron density optimization if `calOption(1)` is true (default: 1).
        - `0`: No minimization (`RhoOptimizers.NoMinimization`).
        - `1`: Sqrt Newton conserving N (`RhoOptN.ProjNormRhoMinimization` with Truncated Newton direction).
        - `2`: CG conserving N (`RhoOptN.ProjNormRhoMinimization` with CG direction).
        - `3`: BFGS conserving N (`RhoOptN.ProjNormRhoMinimization` with BFGS direction).
        - `4`: Sqrt Truncated Newton (`RhoOptSTN.SqrtNewtonMinimization`).
        - `5`: Sqrt Conjugate Gradient (`RhoOptSCG.SqrtGradientMinimization` with `cgmin=.TRUE.`).
        - `6`: Hybrid SCG then STN (`RhoOptSCG.SqrtGradientMinimization` then `RhoOptSTN.SqrtNewtonMinimization`).
        - `7`: Log Truncated Newton (`RhoOptLOG.LogTruncatedNewton`).
    - **`mdType :: INTEGER`**: Selects the type of molecular dynamics simulation (default: -1, meaning no MD).
        - `1`: NVT ensemble (Nose-Hoover thermostat, via `NVT.RunNVT`).
        - `2`: NPT ensemble (Nose-Hoover thermostat + Parrinello-Rahman like barostat, via `NPT.NPTFullCell`).
        - `3`: NVE ensemble (microcanonical, via `NVE.RunNVE`).

- **Main Dispatch Subroutines:**
    - **`SUBROUTINE DoMolecularDynamics`**:
      - **Description:** If `mdType` is set to an MD ensemble type (1, 2, or 3), this routine is called by the main program (`PROFESS.f90`). It further dispatches to the specific MD routine (`RunNVT`, `NPTFullCell`, or `RunNVE`), passing `OptimizeRho` as the callback function for electron density optimization within each MD step. It also handles creating the MD output directory if it doesn't exist.
    - **`SUBROUTINE OptimizeCell`**:
      - **Description:** This routine is called by the main program if `mdType == -1` (i.e., for static calculations).
        - If `calOption(3)` (relax cell) is `.TRUE.`, it calls `SteepestDecent` (from `CellOptimizers`), providing `OptimizeIons` as the callback for optimizing ionic positions within each cell optimization step. The `cellRelax` parameter guides the nature of cell deformation.
        - If `calOption(3)` is `.FALSE.`, it calls `NoOptimization` (from `CellOptimizers`), again passing `OptimizeIons`. This effectively means only ion optimization (if enabled by `calOption(2)`) or a single point calculation is performed. Stress is calculated if `calOption(5)` is true.
    - **`SUBROUTINE OptimizeIons`**:
      - **Description:** This routine is called by `OptimizeCell` (or its delegates). It handles the optimization of ionic positions.
        - If `calOption(2)` (minimize ions) is `.TRUE.`, it selects the ion optimization algorithm based on `ionMethod` and calls the corresponding routine (e.g., `QuickMinOptimization`, `ConjugateGradient`, `IonLBFGS`), passing `OptimizeRho` as the callback for electron density optimization at each new ionic configuration.
        - If `calOption(2)` is `.FALSE.`, it calls `NoOptimization` (from `IonOptimizers`), passing `OptimizeRho`. Forces are calculated if `calOption(4)` is true.
    - **`SUBROUTINE OptimizeRho`**:
      - **Description:** This routine is called by ion/cell optimizers or MD routines. It handles the optimization of the electron density.
        - If `calOption(1)` (minimize density) is `.TRUE.`, it selects the electron density optimization algorithm based on `rhoMethod` and calls the corresponding routine (e.g., `ProjNormRhoMinimization`, `SqrtNewtonMinimization`).
        - If `calOption(1)` is `.FALSE.`, it calls `NoMinimization` (from `RhoOptimizers`) to just calculate the energy for the current fixed density.

# Important Variables/Constants

The module's behavior is dictated by the module-level parameters listed above (`calOption`, `cellRelax`, `ionMethod`, `rhoMethod`, `mdType`). These are typically set during input file processing (e.g., by `ReadInputFile` and stored in `Initializer` or directly here).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

This module acts as a central dispatcher and has numerous dependencies on modules that implement specific algorithms:
- **Task Selection & Control:**
    - `CONSTANTS`: For `DP`.
    - `OUTPUT`, `OutputFiles`: For `errorUnit`, `WrtOut`, `outputUnit`.
    - `SYS`: Provides the data arrays (`rhoR`, `forceIon`, `stress`, `energy`, `frozenIon`, `grids`) that are passed down to the actual worker routines.
    - `MPI_Functions`: For `rankGlobal` (used in `DoMolecularDynamics` for directory creation).
- **Molecular Dynamics:**
    - `MolecularDynamics`: For `md_output_path`.
    - `NVT`: For `RunNVT`.
    - `NPT`: For `NPTFullCell`.
    - `NVE`: For `RunNVE`.
- **Cell Optimization:**
    - `CellOptimizers`: For `NoOptimization` (cell level) and `SteepestDecent`.
- **Ion Optimization:**
    - `IonOptimizers`: For `NoOptimization` (ion level).
    - `IonOptQui`: For `QuickMinOptimization`.
    - `IonOptCG`: For `GradientOptimization` (legacy CG).
    - `IonOptCG2`: For `ConjugateGradient` (revised CG).
    - `IonOptBFGS`: For `IonLBFGS`.
- **Electron Density Optimization:**
    - `RhoOptimizers`: For `NoMinimization` (rho level) and various parameters.
    - `RhoOptSCG`: For `SqrtGradientMinimization`.
    - `RhoOptN`: For `ProjNormRhoMinimization` (which itself uses `RhoDirCG`, `RhoDirBFGS`, `RhoDirNew`).
    - `RhoOptSTN`: For `SqrtNewtonMinimization`.
    - `RhoOptLOG`: For `LogTruncatedNewton`.
- **Timing:**
    - `Timer`: For `StartClock`, `StopClock`.

The main program (`PROFESS.f90`) calls either `DoMolecularDynamics` or `OptimizeCell` based on `mdType`. These routines then cascade calls to `OptimizeIons` and `OptimizeRho` as needed, which in turn select the specific algorithm based on `ionMethod` and `rhoMethod`. This creates a hierarchical control flow for the various computational tasks.
