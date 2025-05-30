# Overview

The `NearestDistance` module (defined in the file `Neighbor.f90`) provides functionality to determine the shortest distance between any two atoms in the simulation cell, taking into account periodic boundary conditions. This check is crucial in simulations like Molecular Dynamics or geometry optimization to prevent atoms from overlapping or becoming unphysically close, which could lead to numerical instability or incorrect physical behavior.

The module offers two distinct methods for this check:
1.  **`NNDist_nobin`**: A direct, brute-force method that calculates the distance between all unique pairs of atoms, considering images in the 26 neighboring periodic cells. This is suitable for smaller systems.
2.  **`NNDist_bin`**: A more computationally efficient method for larger systems that employs a binning (or cell list) technique. The simulation cell is divided into smaller bins, and atoms are sorted into these bins. Distance checks are then primarily performed between atoms within the same bin or in immediately adjacent bins, significantly reducing the number of pairs to check.

The main entry point is `CheckNearestDistanceAtoms`, which selects between these two methods based on the `checkNNDist_bin` flag.

# Key Components

- **`MODULE NearestDistance`**: The main container module.

- **Module-Level Parameters & Variables:**
    - `checkNNDist_bin :: LOGICAL`: A flag that controls the method used for checking the nearest distance. If `.TRUE.` (default), the binning method (`NNDist_bin`) is used. If `.FALSE.`, the direct pairwise check (`NNDist_nobin`) is used.
    - `nnDist :: REAL(kind=DP)`: The minimum allowed distance between any two atoms, in Bohr units (default: 2.0 Bohr). If the actual minimum distance found is less than `nnDist`, the check functions signal an error condition (return -1).

- **`FUNCTION CheckNearestDistanceAtoms(minDistance) RESULT(status)`**:
  - **Description:** This is the primary public interface for the nearest distance check. It calls either `NNDist_bin` or `NNDist_nobin` based on the `checkNNDist_bin` flag.
  - **Arguments:**
    - `minDistance :: REAL(KIND=DP), INTENT(OUT)`: The actual minimum distance found between any pair of atoms.
  - **Return Value:**
    - `status :: INTEGER`: Returns `1` if the minimum distance is greater than or equal to `nnDist` (i.e., no atoms are too close). Returns `-1` if the minimum distance is less than `nnDist`.

- **`FUNCTION NNDist_bin(minDistance) RESULT(status)`**:
  - **Description:** Implements the nearest neighbor distance check using a binning algorithm.
    1.  Skips the check if `nnDist <= 0` or if there's only one atom.
    2.  Defines bin dimensions based on cell vectors and a target `intervalx` (default 20 Bohr).
    3.  Assigns each atom to a bin based on its fractional coordinates.
    4.  Creates a pointer structure (`ptr`, `reg`) to efficiently access atoms within each bin. Atoms are effectively sorted into `sort(3,nions)` array according to their bin.
    5.  Iterates through each bin (`ia,ib,ic`). For each atom in the current bin, it checks distances against other atoms in the same bin and atoms in the 26 neighboring bins (including periodic images by adjusting coordinates if a neighbor bin wraps around the cell).
    6.  Updates `minDistance` if a shorter distance is found.
    7.  Compares the final `minDistance` with `nnDist` to set the return `status`.
  - **Arguments & Return Value:** Same as `CheckNearestDistanceAtoms`.

- **`FUNCTION NNDist_nobin(minDistance) RESULT(status)`**:
  - **Description:** Implements the nearest neighbor distance check by direct pairwise summation.
    1.  Skips the check if `nnDist <= 0` or if there's only one atom.
    2.  Iterates through all unique pairs of atoms (`ia`, `ib` where `ib > ia`).
    3.  For each pair, it calculates the distance between atom `ia` and atom `ib`, and also between atom `ia` and the 26 periodic images of atom `ib` (by adding integer combinations of -1, 0, or 1 times the lattice vectors to `ib`'s coordinates).
    4.  Updates `minDistance` if a shorter distance is found.
    5.  Compares the final `minDistance` with `nnDist` to set the return `status`.
  - **Arguments & Return Value:** Same as `CheckNearestDistanceAtoms`.

# Important Variables/Constants

- **`nnDist`**: The critical threshold distance. If any two atoms are closer than this, it's considered an error or a problematic configuration.
- **`checkNNDist_bin`**: Controls which algorithm is used. The binning algorithm is generally much faster for systems with many atoms.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision) and `PI`.
- **`TIMER`**: For `TimerStart`, `TimerStop`, `stopWatch` type (used for profiling the distance check routines).
- **`CellInfo`**: For the `cell` derived type, which provides access to `cell%ionTable` (ion coordinates) and `cell%cellReal` (real-space lattice vectors).
- **`MATHFUNCTIONS`**: For `Norm` (to calculate vector magnitudes/distances) and `Cross` (used in `NNDist_bin` to determine binning intervals in non-orthogonal cells).
- **`OUTPUT`, `OutputFiles`**: For `WrtOut` (writing messages), `outputUnit`, and the `Title` subroutine (though `Title` is often provided by `MPI_Functions` if it's the same logging utility).
- **`MPI_Functions`**: The `Title` subroutine is imported, potentially for logging. No other direct MPI calls are present in these functions, suggesting they run on each processor with global cell information.

This module is typically called during simulations where atomic positions change (e.g., MD, geometry optimization) to ensure the physical validity of the atomic configuration. If `CheckNearestDistanceAtoms` returns -1, the calling routine might stop the simulation or revert to a previous, valid configuration.
