# Overview

The `ReadIonFile` module is responsible for parsing input files that define the system's geometry, the pseudopotentials for each atomic species, and optionally, an initial electron density. It acts as a bridge between external file formats (often CASTEP-like) and the internal data structures of the PROFESS code, primarily populating the `cell` object (from `CellInfo`) and related arrays like `frozenIon` (from `SYS`).

The module handles:
- Reading cell parameters (lattice vectors or dimensions a,b,c and angles alpha,beta,gamma).
- Reading ion types and their positions (in fractional, Cartesian, or Bohr units).
- Assigning pseudopotential files, atomic core density files (for density decomposition), and atomic density files (for initial guess) to each species.
- Setting up constraints for frozen ions.
- Reading a pre-computed electron density from a file, with capabilities for repeating the unit cell density and interpolating if grid sizes differ.
- Processing pseudopotential files, including unit conversions and setting up data for spline interpolation (used by `LocalPseudoPot`).

# Key Components

- **`MODULE ReadIonFile`**: The main container module.
  - `inputUnit :: INTEGER, PARAMETER`: Default Fortran unit number for reading input files (value: 5).
  - `pseudoFile(:) :: CHARACTER(LEN=systemNameLen + 4), ALLOCATABLE`: A module-level allocatable array used as temporary storage for pseudopotential filenames. It's populated by `ReadGeometry` and then consumed by `ReadPseudo`.

- **`SUBROUTINE ReadGeometry(geometryFile, defaultPseudo)`**:
  - **Description:** Reads the main geometry input file (typically with a `.ion` extension or similar). This file contains information about the simulation cell and the atoms within it.
  - **Parses Blocks:**
    - `%BLOCK LATTICE_CART`: Reads 3x3 real-space lattice vectors (assumed Angstroms, converted to Bohr).
    - `%BLOCK LATTICE_ABC`: Reads a, b, c, alpha, beta, gamma and constructs lattice vectors.
    - `%BLOCK POSITIONS_FRAC` or `POSITIONS_CART` or `POSITIONS_BOHR`: Reads ion species names and their coordinates. Converts Cartesian or Bohr coordinates to fractional if necessary. Sorts ions by type into `cell%ionTable` and populates `cell%elementTable`. Initializes `outputIonTable` for ordered output.
    - `%BLOCK SPECIES_POT`: Assigns pseudopotential filenames to each unique atomic species found. Stores these in the module-level `pseudoFile` array.
    - `%BLOCK SPECIES_CORE`: Assigns atomic core density filenames to species (for `KEDF_DenDec.atomicCoreFile`).
    - `%BLOCK SPECIES_RHOA`: Assigns atomic density filenames to species (for `AtomicDensity.atomicDensityFile`).
    - `%BLOCK ION_OPTIMIZATION`: Reads flags (0 or 1) for each ion and each Cartesian direction to determine if that degree of freedom is frozen (populates `SYS.frozenIon`).
  - **Actions:** Populates `cell%cellReal`, `cell%ionTable`, `cell%elementTable`, `SYS.frozenIon`, `pseudoFile`, `AtomicDensity.atomicDensityFile`, `KEDF_DenDec.atomicCoreFile`. Calls `RefreshLattice`.
  - **Contains Internal Function:**
    - `FUNCTION Locate(x, a) RESULT(index)`: Finds the first index of integer `x` in integer array `a`.

- **`SUBROUTINE ReadPseudo`**:
  - **Description:** Reads and processes the pseudopotential files for each unique atomic species identified in `ReadGeometry`. The filenames are taken from the module-level `pseudoFile` array.
  - **For each pseudopotential file (typically `.recpot`):**
    1.  Opens the file.
    2.  Determines type (currently only `.recpot` seems fully supported).
    3.  Skips comment lines (ending with "END COMMENT") and version information.
    4.  Reads `psp%maxG` (max G-vector magnitude, converted from A^-1 to Bohr^-1) and the number of data points `psp%numPoints`.
    5.  Reads the tabulated pseudopotential values `V(q)`.
    6.  Converts `V(q)` units from `eV*AU^3/A^3` (likely a typo in comments, probably `eV*Angstrom^3`) to atomic units (Hartree*Bohr^3).
    7.  Calculates the ionic charge `Z_eff` from the behavior of `V(q)` near `q=0` (using `V(q) ~ -4*pi*Z_eff/q^2 + V(q=0)`).
    8.  Stores `Z_eff` in `psp%charge` and updates `cell%elementTable(i)%charge` and `cell%elementTable(i)%chargeTot`.
    9.  Prepares data for spline interpolation:
        - Stores `q` values (as indices) in `psp%t`.
        - Stores `V(q) + 4*pi*Z_eff/q^2` in `psp%vqS` (this quantity is finite at q=0).
        - Calls `spline_cubic_set` (from `MathSplines`) to compute and store second derivatives of `psp%vqS` into `psp%potDD`.
    10. Assigns the populated `psp` pointer to `cell%elementTable(i)%psp`.
  - **Side Effects:** Populates `cell%elementTable(i)%psp` for all `i`, calculates and stores `cell%numEle`. Deallocates `pseudoFile`.

- **`SUBROUTINE ReadDensity(densityFile, xRepeat, yRepeat, zRepeat)`**:
  - **Description:** Reads a 3D real-space electron density from a specified `densityFile` into the global `SYS.rhoR` array. The input file is expected to be a direct access, formatted file with a header containing dimensions (`xLen, yLen, zLen`) and number of spin channels (`fileSpin`).
  - **Functionality:**
    - If the dimensions in the file (optionally repeated by `xRepeat, yRepeat, zRepeat`) match the current simulation grid dimensions, the data is copied directly.
    - If dimensions differ, trilinear interpolation is performed to map the density from the file grid onto the simulation grid.
    - After reading/interpolating, the total electron density is normalized to match the system's expected total number of electrons (`numEle` derived from ionic charges).
  - **Arguments:**
    - `densityFile :: CHARACTER(LEN=*), INTENT(IN)`: Path to the density file.
    - `xRepeat, yRepeat, zRepeat :: INTEGER, INTENT(IN)`: Number of times to repeat the unit cell density from the file in each direction.

# Important Variables/Constants

- **`inputUnit :: INTEGER, PARAMETER = 5`**: Default Fortran unit for reading these files.
- **`pseudoFile(:)`**: Temporary storage for pseudopotential filenames between `ReadGeometry` and `ReadPseudo`.
- **File Format Assumptions:** The routines assume specific block structures (e.g., `%BLOCK LATTICE_CART`) and data orders within the input files, often similar to CASTEP formats. Pseudopotential files also have an assumed structure.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`, `systemNameLen`, `PI`, `bohr`, `hartreeToeV`, `AtomicMass`.
- **`Output`, `OutputFiles`**: For `WrtOut`, `QUIT`, `errorUnit`, `outputUnit`, `Error`, and `outputIonTable`.
- **`CellInfo`**: For the `cell` derived type (which is heavily populated by this module), `ion` type, `pseudoPot` type, and `RefreshLattice`.
- **`Sys`**: For `frozenIon` array (populated by `ReadGeometry`) and `rhoR` array (populated by `ReadDensity`).
- **`MathFunctions`**: For `Volume`, `Cross`, `Vecmul`, `Inverse`.
- **`MathSplines`**: For `spline_cubic_set` (used in `ReadPseudo` to prepare for pseudopotential interpolation).
- **`AtomicDensity`**: The `atomicDensityFile` array (module variable in `AtomicDensity`) is populated by `ReadGeometry`.
- **`KEDF_DenDec`**: The `atomicCoreFile` array (module variable in `KEDF_DenDec`) is populated by `ReadGeometry`.
- **`MPI_Functions`**: For `ReduceRealLevel1` and `mpiErr` (used in `ReadDensity` for parallel normalization).

This module is called early in the simulation setup process by `Initializer.InitPROFESS` to define the physical system based on user-provided files. The data structures it populates are fundamental to almost all subsequent calculations.
