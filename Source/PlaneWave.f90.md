# Overview

The `PlaneWave` module is responsible for managing fundamental quantities related to the plane wave (reciprocal space, G-space) representation used in electronic structure codes. Its primary tasks include:
1.  Generating the set of reciprocal lattice vectors (G-vectors or q-vectors) based on the simulation cell's reciprocal lattice and a specified kinetic energy cutoff.
2.  Storing the magnitudes (norms) of these G-vectors in `qTable`.
3.  Storing the Cartesian components of these G-vectors in `qVectors`.
4.  Creating a boolean mask (`qMask`) that identifies which G-vectors are within the energy cutoff and are part of the non-redundant half-sphere used in real-to-complex Fast Fourier Transforms (FFTs).
5.  Providing a utility function `CCStructureFactor` to calculate the complex conjugate of the ionic structure factor for a given set of ions and a G-vector.
6.  Offering a helper function `WrapQ` to adapt a G-space table to a different (typically smaller) grid size, although its implementation appears to be a simple truncation/mirroring.

# Key Components

- **`MODULE PlaneWave`**: The main container for plane wave related data and routines.

- **Module-Level Variables (Allocatable):**
    - **`qTable(:,:,:) :: REAL(KIND=DP)`**: Stores the norm (magnitude, `|G|`) of each reciprocal lattice vector `G` corresponding to a point on the G-space FFT grid.
    - **`qVectors(:,:,:,:) :: REAL(KIND=DP)`**: Stores the Cartesian components `(Gx, Gy, Gz)` of each G-vector. The dimensions are typically `(k1G, k2G, k3G, 3)`, where `k1G, k2G, k3G` are the dimensions of the G-space grid for the current MPI process.
    - **`qMask(:,:,:) :: LOGICAL`**: A boolean mask corresponding to the G-space grid. An element `qMask(i,j,k)` is `.TRUE.` if the G-vector `qVectors(i,j,k,:)` is within the specified kinetic `energyCutoff` AND it belongs to the non-redundant set of G-vectors (typically a half-sphere, due to `V(G) = V(-G)^*` for real potentials/densities). The `G=0` point is generally excluded by this mask (i.e., `qMask(1,1,1)` is `.FALSE.`).

- **Module-Level Parameter:**
    - `energyCutoff :: REAL(KIND=DP)`: The kinetic energy cutoff (e.g., in eV or Hartree) that defines the maximum extent of G-vectors included in calculations. `E_cut = (hbar^2 / 2m) * G_max^2`. Default is -1, indicating it must be set.

- **`FUNCTION CCStructureFactor(singleTypeIonTable, cellReal, qPoint) RESULT(structure_factor_conj)`**:
  - **Description:** Calculates the complex conjugate of the structure factor for a specific ion type. The formula is `S_type(G)^* = SUM_j { exp(i * G . r_j) }`, where the sum is over all ions `j` of the given type, and `r_j` are their real-space positions.
  - **Arguments:**
    - `singleTypeIonTable :: TYPE(ion), DIMENSION(:), INTENT(IN)`: Array of ion data for a single element type.
    - `cellReal :: REAL(KIND=DP), DIMENSION(3,3), INTENT(IN)`: Real-space lattice vectors.
    - `qPoint :: REAL(KIND=DP), DIMENSION(3), INTENT(IN)`: The G-vector (Cartesian coordinates) at which to calculate the structure factor.
  - **Return Value:** `structure_factor_conj :: COMPLEX(KIND=DP)`.

- **`SUBROUTINE FillQTable(cellRecip)`**:
  - **Description:** Populates the module-level arrays `qVectors`, `qTable`, and `qMask`. It iterates through all points on the G-space FFT grid (as defined by `k1G, k2G, k3G, k3Goff` from `CellInfo` and global dimensions `m1G, m2G, m3G`). For each grid point, it constructs the corresponding G-vector using integer multiples of reciprocal lattice basis vectors (`cellRecip`), calculates its norm (stored in `qTable`), and stores its Cartesian components (in `qVectors`). The `qMask` is set based on whether `|G|^2` is within the `energyCutoff` and based on symmetry considerations for real-to-complex FFTs (selecting roughly half of the G-space points).
  - **Arguments:**
    - `cellRecip :: REAL(kind=DP), DIMENSION(:,:), INTENT(IN)`: The reciprocal space lattice vectors.

- **`FUNCTION WrapQ(table, x, y, z) RESULT(wrapped_table)`**:
  - **Description:** A utility function intended to map or interpolate a G-space table (like `qTable`) onto a grid of a different size (`x,y,z`), typically if the new grid is smaller or the same. The current implementation seems to perform a direct copy for the overlapping region and then uses mirroring for points outside the original table's bounds but within the new smaller bounds. It does not perform sophisticated interpolation.
  - **Arguments:**
    - `table :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: The input G-space table.
    - `x, y, z :: INTEGER, INTENT(IN)`: Dimensions of the target wrapped table.
  - **Return Value:** `wrapped_table :: REAL(KIND=DP), DIMENSION(x,y,z)`.

# Important Variables/Constants

- **`energyCutoff`**: Defines the sphere in reciprocal space containing all G-vectors to be considered.
- **`qTable`**: Stores `|G|`. Essential for functions that depend on the magnitude of G, like radial pseudopotentials or some kernel functions.
- **`qVectors`**: Stores `(Gx, Gy, Gz)`. Used for dot products (e.g., `G.r` in structure factors) and terms involving components of G.
- **`qMask`**: Crucial for summations in reciprocal space, ensuring only relevant and non-redundant G-vectors are included.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision), `PI`, and `IMAG` (imaginary unit).
- **`OutputFiles`**: For `outputUnit`.
- **`FOURIER`**: For `offset` (an integer related to FFT grid parity, likely `MOD(global_dimX_real, 2)`), which is used in `FillQTable` to correctly set up `qMask` for different FFT symmetries.
- **`MathFunctions`**: For `Vecmul` (matrix-vector multiplication, used to get Cartesian G-vectors from fractional multiples of reciprocal lattice vectors) and `Norm` (to calculate `|G|`).
- **`CellInfo`**: For `ion` type definition (used in `CCStructureFactor`), and for various G-space and real-space FFT grid dimensions (`k1G, k2G, k3G, k3Goff, m1G, m2G, m3G`) used in `FillQTable`.
- **`MPI_Functions`**: The `TITLE` subroutine is imported but not directly used in the provided code snippets.

The `PlaneWave` module is fundamental for any calculations performed in reciprocal space. `FillQTable` is typically called early in the setup of a calculation (e.g., by `Sys` or `Initializer` modules, after cell parameters are known and FFTs are planned) to make `qTable`, `qVectors`, and `qMask` available to all other modules that need them (e.g., KEDF modules, Hartree potential calculation, pseudopotential evaluation). `CCStructureFactor` is used by modules like `IonElectron` and `KEDF_DenDec` when constructing potentials or densities from atomic contributions.
