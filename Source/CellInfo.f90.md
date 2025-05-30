# Overview

The `CellInfo` module serves as a central repository for crucial information about the simulation cell, the ions it contains, and the pseudopotentials used. It defines several data types (`pseudoPot`, `element`, `ion`, `cellStruct`) to organize this information. Additionally, it stores parameters related to Fast Fourier Transform (FFT) grids in both real and reciprocal space. The module provides subroutines to initialize these FFT dimensions and to update lattice-derived quantities (like reciprocal vectors, volume, and lattice parameters) whenever the simulation cell geometry changes.

# Key Components

- **`MODULE CellInfo`**: The main container for all data types, module variables, and subroutines.

- **Data Types:**
    - **`TYPE pseudoPot`**: Represents a pseudopotential.
        - `numPoints :: INTEGER`: Total number of points in the pseudopotential data.
        - `charge :: INTEGER`: Ionic charge associated with the pseudopotential (Z_ion - Z_core).
        - `maxG :: REAL(KIND=DP)`: Maximum wavevector (G) if in reciprocal space (.recpot).
        - `maxR :: REAL(KIND=DP)`: Maximum radius (r) if in real space (.realpot).
        - `potValues :: REAL(KIND=DP), DIMENSION(:), POINTER`: Pointer to array of pseudopotential values.
        - `potDD :: REAL(KIND=DP), DIMENSION(:), POINTER`: Pointer to array of second derivatives (for spline interpolation).
        - `vqS :: REAL(KIND=DP), DIMENSION(:), POINTER`: Pointer to array of `v(q) + 4*pi*Charge/q**2` (for spline interpolation).
        - `t :: REAL(KIND=DP), DIMENSION(:), POINTER`: Pointer to array of knots for spline interpolation.
        - `pspFileName :: CHARACTER(systemNameLen+7)`: Filename of the pseudopotential.
        - `type :: INTEGER`: Type of pseudopotential (1 for `.recpot`, 2 for `.realpot`).

    - **`TYPE element`**: Represents an atomic element type.
        - `elementName :: CHARACTER(LEN=2)`: Chemical symbol (e.g., "Si").
        - `charge :: REAL(KIND=DP)`: Effective charge per atom of this type.
        - `chargeTot :: REAL(KIND=DP)`: Total electronic charge for all atoms of this type.
        - `mass :: REAL(KIND=DP)`: Atomic mass.
        - `psp :: TYPE(pseudoPot), POINTER`: Pointer to the `pseudoPot` for this element.
        - `firstIonID :: INTEGER`: Index of the first ion of this type in the `ionTable`.

    - **`TYPE ion`**: Represents an individual ion.
        - `elementID :: INTEGER`: Identifier linking to the `elementTable`.
        - `coord :: REAL(kind=DP), DIMENSION(3)`: Fractional coordinates of the ion.

    - **`TYPE cellStruct`**: Represents the simulation cell and its contents.
        - `cellReal :: REAL(KIND=DP), DIMENSION(3,3)`: Real-space lattice vectors (in Bohr).
        - `cellRecip :: REAL(KIND=DP), DIMENSION(3,3)`: Reciprocal-space lattice vectors.
        - `lengthA, lengthB, lengthC :: REAL(KIND=DP)`: Magnitudes of the real-space lattice vectors.
        - `angleAlpha, angleBeta, angleGamma :: REAL(KIND=DP)`: Angles between lattice vectors (in degrees).
        - `vol :: REAL(KIND=DP)`: Cell volume (in Bohr^3).
        - `dv :: REAL(KIND=DP)`: Real-space grid element volume (dV = vol / m123G).
        - `numIonType :: INTEGER`: Number of distinct element types.
        - `numIon :: INTEGER`: Total number of ions in the cell.
        - `numEle :: REAL(KIND=DP)`: Total number of electrons.
        - `ionTable :: TYPE(ion), DIMENSION(:), POINTER`: Pointer to an array of `ion` types.
        - `elementTable :: TYPE(element), DIMENSION(:), POINTER`: Pointer to an array of `element` types.
        - `totNumEle :: REAL(kind=DP)`: (Seems redundant with `numEle`) Total number of electrons.

- **Subroutines:**
    - **`SUBROUTINE SetupCellFFTdims`**:
        - **Description:** Initializes the local and global FFT grid dimensions for both real space (`n1G`, `n2G`, `n3G`, `m1G`, `m2G`, `m3G`) and reciprocal space (`k1G`, `k2G`, `k3G`). Also calculates total grid points (`n123G`, `m123G`, `k123G`) and local offsets (`n3Goff`, `k3Goff`) for parallel execution.
    - **`SUBROUTINE RefreshLattice(newCell)`**:
        - **Description:** Updates various cell-derived quantities when the cell geometry changes. Given a new set of real-space lattice vectors (`newCell`), it recalculates `cell%cellReal`, `cell%cellRecip`, `cell%vol`, `cell%dv`, and the lattice parameters `lengthA, lengthB, lengthC, angleAlpha, angleBeta, angleGamma`.
        - **Arguments:**
            - `newCell :: REAL(KIND=DP), DIMENSION(3,3), INTENT(IN)`: The new real-space cell matrix.

# Important Variables/Constants (Module-Level)

**FFT Grid Parameters (Real Space):**
- `n1G, n2G, n3G :: INTEGER`: Local FFT dimensions on the current processor.
- `n123G :: INTEGER`: Total number of local real-space grid points (`n1G*n2G*n3G`).
- `m1G, m2G, m3G :: INTEGER`: Global FFT dimensions for the entire cell.
- `m123G :: INTEGER`: Total number of global real-space grid points (`m1G*m2G*m3G`).
- `n3Goff :: INTEGER`: Offset for the local grid in the z-dimension (for slab decomposition in parallel FFT).

**FFT Grid Parameters (Reciprocal Space):**
- `k1G, k2G, k3G :: INTEGER`: Local G-space (reciprocal space) FFT dimensions.
- `k123G :: INTEGER`: Total number of local G-space grid points (`k1G*k2G*k3G`).
- `k3Goff :: INTEGER`: Offset for the local G-space grid in the z-dimension.

**Other Module Variables:**
- `numSpin :: INTEGER`: Number of electron spin channels (default is 1 for spin-unpolarized, 2 for spin-polarized).
- `cell :: TYPE(cellStruct), SAVE`: The primary variable instance of `cellStruct`, holding all geometric and compositional information about the system. `SAVE` attribute ensures it persists throughout the program execution.
- `usingGridPack :: LOGICAL`: Flag to indicate if grid packing for multigrid methods is active (default `.FALSE.`).
- `storeRealIonPot :: LOGICAL`: Flag indicating whether the real-space ionic potential should be stored (default `.TRUE.`, often required for multigrid).
- `kinetic :: INTEGER`: Default selector for the Kinetic Energy Density Functional (KEDF). Initialized to 5, which typically corresponds to the WGC KEDF in this codebase.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter), `systemNameLen` (length for character strings like filenames), and `PI`.
- **`OUTPUTFILES`**: For `outputUnit` (Fortran I/O unit number for standard output).
- **`FOURIER`**: For `GetFFTDims` and `GetFFTComplexDims` subroutines, which are used in `SetupCellFFTdims` to determine the FFT grid layout based on the parallel configuration.
- **`MathFunctions`**: For `Volume` (calculates cell volume from lattice vectors) and `Inverse` (calculates the inverse of a matrix), used in `RefreshLattice`.
- **`MPI_Functions`**: For `message`, `Error` (error handling routines), and `Title` (routine to print titles/headers, though this might also be from an `Output` module).

The `CellInfo` module is fundamental to the application, providing the data structures and basic geometric information required by almost all other physics and calculation modules. It interacts closely with input routines (which populate its structures) and with any routine that needs to know about cell dimensions, ion positions, or FFT grid details.
