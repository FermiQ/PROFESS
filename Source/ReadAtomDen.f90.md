# Overview

The `AtomicDensity` module (defined in the file `ReadAtomDen.f90`) is responsible for constructing an initial guess for the total electron density of the system. It does this by superimposing pre-calculated atomic electron densities.

The process involves:
1.  Reading 1D atomic charge densities in reciprocal space (`rho_atom(G)`) from specified files for each unique atomic species in the system.
2.  For each point `G` on the 3D reciprocal space grid of the simulation cell:
    a.  Interpolating the 1D `rho_atom(|G|)` for each atomic species to get its value at `|G|`.
    b.  Multiplying this interpolated atomic density by the appropriate ionic structure factor `S(G)` for that species.
    c.  Summing these contributions over all species to get the total `rho_atomic_sum(G)`.
3.  Performing an inverse Fast Fourier Transform (FFT) on `rho_atomic_sum(G)` to obtain the total initial electron density in real space, `rho(r)`.

This initial density is often used as a starting point for self-consistent field (SCF) calculations or direct minimization methods.

# Key Components

- **`MODULE AtomicDensity`**: The main container module.

- **`TYPE atomDen`**:
  - **Description:** A derived data type to store the tabulated 1D atomic density data read from a file.
  - **Fields:**
    - `nG :: INTEGER`: The total number of points in the tabulated 1D atomic density.
    - `maxG :: REAL(KIND=DP)`: The maximum G-vector magnitude for which the 1D atomic density is defined.
    - `rhoG(:) :: REAL(KIND=DP), ALLOCATABLE`: Array of atomic density values `rho_atom(G)` in reciprocal space.
    - `normG(:) :: REAL(KIND=DP), ALLOCATABLE`: Array of corresponding G-vector magnitudes.

- **Module-Level Variables:**
    - `atomicDensityFile(:) :: CHARACTER(LEN=100), ALLOCATABLE`: An array storing the filenames for the 1D atomic density data for each atomic species. This array is expected to be allocated and populated externally (e.g., by `ReadInputFile` based on user input).
    - `atomTypeNum :: INTEGER`: The number of different atomic types for which atomic density files are provided.

- **`SUBROUTINE GenerateAtomicDensity(rho, dimX, dimY, dimZ, nspin)`**:
  - **Description:** This is the main public subroutine. It orchestrates the generation of the total real-space atomic density.
    1.  Determines `atomTypeNum` from `SIZE(atomicDensityFile)`.
    2.  Allocates an array `aden(atomTypeNum)` of type `atomDen`.
    3.  Calls the internal `ReadAtomicDensity` subroutine to populate the `aden` array.
    4.  For each spin channel (`is = 1, nspin`), it calls the internal `ADenReal` subroutine to compute the real-space atomic density and store it in the corresponding spin channel of the output array `rho`.
  - **Arguments:**
    - `rho :: REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(INOUT)`: The output 3D real-space electron density array (with spin dimension).
    - `dimX, dimY, dimZ, nspin :: INTEGER, INTENT(IN)`: Dimensions of the `rho` array.

- **`SUBROUTINE ReadAtomicDensity`** (Internal to `GenerateAtomicDensity`):
  - **Description:** Reads the 1D atomic density data (`|G|` and `rho_atom(G)`) from the files specified in `atomicDensityFile(:)` for each atom type. The data is stored in the `aden` array of `atomDen` type. It converts the read-in `|G|` values (assumed to be in Angstrom^-1) to Bohr^-1 by multiplying by `bohr`.
  - **Output:** Populates the `aden` array (module variable accessible within `GenerateAtomicDensity`).

- **`FUNCTION ADenLookup(aden_type, qNorm) RESULT(interpolated_rhoG)`** (Private):
  - **Description:** Interpolates the 1D atomic density `aden_type%rhoG` (which is `rho_atom(|G|)`) to find its value at a specific `qNorm` (magnitude of a G-vector). It uses a 4-point Lagrange-like interpolation formula. Special handling is provided for `qNorm = 0` (returns `aden_type%rhoG(1)`) and `qNorm >= aden_type%maxG` (returns 0.0).
  - **Arguments:**
    - `aden_type :: TYPE(atomDen), INTENT(IN)`: The `atomDen` object for a specific atom type.
    - `qNorm :: REAL(KIND=DP), INTENT(IN)`: The G-vector magnitude at which to interpolate.
  - **Return Value:** `interpolated_rhoG :: REAL(KIND=DP)`.

- **`FUNCTION ADenRecip(ionTable, elementTable, aden) RESULT(total_atomic_rho_recip)`** (Private):
  - **Description:** Constructs the total 3D reciprocal space atomic density `rho_atomic_sum(G)`. It iterates over all G-vectors on the simulation grid. For each `G`:
    1.  Calculates its norm, `|G|`.
    2.  For each atom type, calls `ADenLookup` to get `rho_atom(|G|; type)`.
    3.  Multiplies this by the complex conjugate of the structure factor `CCStructureFactor(type, G)`.
    4.  Sums these contributions over all atom types.
    5.  The final sum is divided by the cell volume.
  - **Arguments:**
    - `ionTable :: TYPE(ion), DIMENSION(:), INTENT(IN)`: System's ion data.
    - `elementTable :: TYPE(element), DIMENSION(:), INTENT(IN)`: System's element data.
    - `aden :: TYPE(atomDen), DIMENSION(:), INTENT(IN)`: Array of loaded atomic densities.
  - **Return Value:** `total_atomic_rho_recip :: COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G)`.

- **`SUBROUTINE ADenReal(rho_out, aden_in)`** (Private):
  - **Description:** Calculates the real-space atomic density `rho_out` by performing an inverse FFT (using `FFT_NEW` via the generic `FFT` interface) on the reciprocal space density obtained from `ADenRecip`.
  - **Arguments:**
    - `rho_out :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(INOUT)`: Output real-space density.
    - `aden_in :: TYPE(atomDen), DIMENSION(:), INTENT(IN)`: Array of loaded atomic densities (passed to `ADenRecip`).

# Important Variables/Constants

- **`atomicDensityFile(:)`**: Array of character strings holding the paths to the files containing 1D atomic density data. This is the primary input defining where to get atomic data.
- **`aden(:)` (within `GenerateAtomicDensity`)**: An array of `atomDen` type, storing the loaded 1D atomic densities for each species.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision), `PI`, and `bohr` (for unit conversion).
- **`OutputFiles`**: For `outputUnit` (Fortran I/O unit for messages).
- **`PlaneWave`**: For `qMask` (imported but not directly used in the core logic shown), `qVectors` (G-vector components), and `CCStructureFactor` (ionic structure factor).
- **`CellInfo`**: For `cell` derived type (providing `ionTable`, `elementTable`, `cell%cellReal`, `cell%vol`), and G-space grid dimensions (`k1G, k2G, k3G`).
- **`MathFunctions`**: For `Volume` (cell volume) and `Norm` (magnitude of G-vectors).
- **`Fourier_NEW`**: For `FFT_NEW` and `FFT_STD_STATE`, which are used by `ADenReal` via the generic `FFT` interface provided by the `Fourier` shim module (though `Fourier` itself is not directly `USE`d here, `FFT_NEW` implies its usage context).
- **`MPI_Functions`**: For `Title`, `StartClock`, `StopClock` (logging and timing utilities).

This module is typically called during the initialization phase of a PROFESS calculation (e.g., from the `Initializer` module) if the input options specify that the starting electron density should be constructed from atomic densities. The resulting `rho` array is then used by other parts of the code, such as the SCF cycle or direct minimization routines.
