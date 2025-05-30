# Overview

The `CellOptimizers` module provides subroutines for optimizing the simulation cell's lattice vectors based on the calculated stress tensor. These routines typically work in conjunction with an external ion/electron structure optimizer. After each modification to the cell, this external optimizer is called to relax the ionic positions and electron density within the new cell geometry before the stress is recalculated.

# Key Components

- **`MODULE CellOptimizers`**:
  - The main container for cell optimization routines.
  - Defines module-level parameters:
    - `maxCellStep :: INTEGER`: Maximum number of cell optimization steps allowed (default: 100).
    - `tols :: REAL(KIND=DP)`: Tolerance for the stress tensor components to determine convergence (default: 5.0E-7 a.u.).

- **`SUBROUTINE NoOptimization(Optimizer, rhoR, energy, stress, calcStress)`**:
  - **Description:** This routine does not perform any cell optimization. It simply calls the provided `Optimizer` subroutine (which is expected to handle ion and electron relaxation) and then, if `calcStress` is true, it calls `CalculateStress` to compute the stress tensor for the current (fixed) cell.
  - **Arguments:**
    - `Optimizer :: EXTERNAL`: An external subroutine that optimizes ion positions and electron density for the given cell.
    - `rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: The electron density.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: Table of energy components.
    - `stress :: REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT)`: The calculated stress tensor.
    - `calcStress :: LOGICAL, INTENT(IN)`: Flag to indicate whether to calculate the stress.

- **`SUBROUTINE SteepestDecent(Optimizer, rhoR, energy, stress, axis)`**:
  - **Description:** Implements a steepest descent-like algorithm to minimize stress by adjusting the cell lattice vectors. It iteratively updates the cell based on the current stress, calls the `Optimizer` to relax ions/electrons, recalculates stress, and checks for convergence. The step size `dt` is adapted based on changes in total energy. It can perform optimization along a specific axis (x, y, or z) or for all cell parameters.
  - **Arguments:**
    - `Optimizer :: EXTERNAL`: External subroutine for ion/electron relaxation.
    - `rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`: The real-space electron density (modified by `RescaleDensity`).
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: Table of energy components; `energy(1)` (total energy) is monitored.
    - `stress :: REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT)`: The calculated stress tensor.
    - `axis :: INTEGER, INTENT(IN)`: Specifies the optimization mode:
        - `1`: Optimize only the x-direction components of stress/cell.
        - `2`: Optimize only the y-direction.
        - `3`: Optimize only the z-direction.
        - `Otherwise`: Full cell optimization.
  - **Internal Subroutine:**
    - **`SUBROUTINE WriteLogFile`**: Writes detailed log information about the current cell optimization step (step number, max stress, total energy, stress tensor, lattice vectors) to the standard output and potentially other log files.

# Important Variables/Constants

- **`maxCellStep :: INTEGER`**: Module-level parameter defining the maximum iterations for cell optimization loops.
- **`tols :: REAL(KIND=DP)`**: Module-level parameter for the convergence threshold of the stress tensor.
- **`Optimizer :: EXTERNAL`**: A placeholder for an externally provided subroutine responsible for relaxing ion positions and electron density.
- **`rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:)`**: The electron density. In `SteepestDecent`, it's `INTENT(INOUT)` because it's rescaled after cell changes.
- **`energy :: REAL(KIND=DP), DIMENSION(:)`**: Array of energy components. `energy(1)` is used in `SteepestDecent` to adapt the step size.
- **`stress :: REAL(KIND=DP), DIMENSION(3,3)`**: The stress tensor.
- **`axis :: INTEGER`**: Input to `SteepestDecent` to control the type of cell optimization (uniaxial or full).
- **`dt :: REAL(KIND=DP)`**: Internal variable in `SteepestDecent` representing the adaptive step size for updating the cell matrix.
- **`newCell :: REAL(KIND=DP), DIMENSION(3,3)`**: Internal variable in `SteepestDecent` to temporarily hold the modified cell lattice vectors.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- **`CellInfo`**:
    - For the `cellStruct` derived data type.
    - For the `cell` module variable (instance of `cellStruct`).
    - For `RefreshLattice` subroutine (to update cell parameters after modification).
- **`RefreshCell`**: For `RefreshCellSetup` subroutine (to update other system properties that depend on cell size/shape after modification).
- **`CalStress`**: For `CalculateStress` subroutine (to compute the stress tensor).
- **`RefreshIons`**: For `RescaleDensity` subroutine (to ensure electron density integrates to the correct number of electrons after cell volume changes).
- **`OutputFiles`**: For `outputUnit` (Fortran I/O unit number).
- **`MathFunctions`**: For `Inverse` function (to calculate the inverse of the cell matrix).
- **`Output`**:
    - For `WrtOut` (general output subroutine).
    - For `PrintStress` (to print the stress tensor).
    - For `PrintGeometry` (to print cell geometry).
    - For `celOutputFreq` (controls frequency of geometry printing).

The `CellOptimizers` module acts as a higher-level manager for cell relaxation. It relies on:
1.  An external `Optimizer` routine to handle the electronic and ionic degrees of freedom for a *given* cell.
2.  The `CalStress` module to calculate the stress tensor for that relaxed state.
3.  Various routines from `CellInfo`, `RefreshCell`, and `RefreshIons` to update the system's description and state when the cell itself is changed.
