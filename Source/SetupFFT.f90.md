# Overview

The `SetupFFT` module is responsible for initializing and configuring the Fast Fourier Transform (FFT) grids and plans for the entire simulation. It determines the optimal dimensions of the FFT grid based on either a specified kinetic energy cutoff (`energyCutoff`) or a direct grid spacing (`gridSpacing`). Once dimensions are determined, it calls routines from the `FOURIER_NEW` module (which interfaces with FFTW3) to prepare the FFT plans.

The module handles both serial and parallel FFT setup. In parallel, it uses a specific MPI communicator (`MPI_GFFT_WORLD`) for planning global FFTs and then retrieves local (per-processor) FFT dimensions and offsets.

# Key Components

- **`MODULE SetupFFT`**: The main container module.

- **Module-Level Parameter:**
    - `dimType :: INTEGER`: Controls the scheme for padding FFT dimensions in `PadFFTdims`. Default is 1.
        - `1`: Dimensions are products of small primes (2, 3, 5, 7).
        - `2`: Same as 1, but ensures at least one factor of 2 (even dimension).
        - `3`: Dimensions are products of odd primes (3, 5, 7).
        - `4`: Dimensions are powers of 2.
        - `5`: Dimensions are powers of 2 multiplied by at most one factor of 3, 5, or 7.

- **`SUBROUTINE InitializeFFT(energyCutoff, gridSpacing, totX, totY, totZ, locZ, locZOff)`**:
  - **Description:** This is the primary routine for setting up the FFTs.
    1.  Determines global FFT dimensions (`totX`, `totY`, `totZ`):
        - If `energyCutoff > 0.0`, it calls `SizeSystem` to calculate dimensions based on the cutoff and cell parameters.
        - Otherwise, it calculates dimensions based on `gridSpacing` and cell lengths.
    2.  Prints the global FFT dimensions.
    3.  Calls `PlanFFT_NEW` (from `FOURIER_NEW`) to create FFT plans. By default, it uses `FFTW3_ONEPLAN_ESTIMATE` for planning strategy. For parallel runs (`__USE_PARALLEL` defined), it passes the `MPI_GFFT_WORLD` communicator.
    4.  If running in parallel:
        - Calls `GetFFTComplexDims` and `GetFFTDims` (from the `FOURIER` shim module, which in turn call `_NEW` versions) to obtain the local z-dimension (`locZ`) and offset (`locZOff`) for the current processor's slab of the FFT grid.
        - Prints these local dimensions.
        - Includes a check to ensure consistency between `locZOff` and `rankGFFT * locZ`.
  - **Arguments:**
    - `energyCutoff :: REAL(KIND=DP), INTENT(IN)`: Kinetic energy cutoff to determine grid size.
    - `gridSpacing :: REAL(KIND=DP), INTENT(IN)`: Alternative way to define grid size if `energyCutoff <= 0`.
    - `totX, totY, totZ :: INTEGER, INTENT(INOUT)`: Global FFT dimensions (output if calculated from `energyCutoff`, input if from `gridSpacing` but potentially modified by padding).
    - `locZ, locZOff :: INTEGER, INTENT(INOUT)`: Local Z-dimension and offset for parallel FFTs (output).

- **`SUBROUTINE SizeSystem(energyCutoff, numProc, numX, numY, numZ)`**:
  - **Description:** Calculates the required global FFT grid dimensions (`numX`, `numY`, `numZ`) based on the `energyCutoff` and the cell's lattice vector lengths. The formula `qMax = SQRT(energyCutoff * 2.0)` determines the maximum reciprocal space vector magnitude, which then translates to minimum real-space grid points `m?Max = qMax * L? / PI`. These values are then passed to `PadFFTdims`.
  - **Arguments:**
    - `energyCutoff :: REAL(KIND=DP), INTENT(IN)`.
    - `numProc :: INTEGER, INTENT(IN)`: Number of processors in the FFT group (used by `PadFFTdims` if `procMultiple` is true).
    - `numX, numY, numZ :: INTEGER, INTENT(OUT)`: Calculated global grid dimensions.
  - **Contains Internal Function:**
    - **`FUNCTION PadFFTdims(R, procMultiple) RESULT(padded_dim)`**:
      - **Description:** Takes a real number `R` (representing a minimum required grid dimension) and returns an integer `padded_dim` that is greater than or equal to `R`. The returned dimension is chosen to be "FFT-friendly" by being composed of products of small prime numbers (2, 3, 5, 7), according to the module-level `dimType` parameter. If `procMultiple` is `.TRUE.`, it also ensures `padded_dim` is a multiple of `numProc` (relevant for parallel FFTs where `numProc` is `sizeGFFT`).
      - **Arguments:**
        - `R :: REAL(KIND=DP), INTENT(IN)`: The minimum dimension required.
        - `procMultiple :: LOGICAL, INTENT(IN)`: If true, ensures the result is a multiple of `numProc`.

# Important Variables/Constants

- **`dimType :: INTEGER`**: Module-level parameter controlling the padding strategy for FFT dimensions to ensure efficient FFT execution.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision) and `PI`.
- **`MPI_Functions`**: For MPI-related variables: `rankGlobal` (current process rank in `MPI_COMM_WORLD`), `sizeGFFT` (number of processors in the FFT group), `MPI_GFFT_WORLD` (the communicator for the FFT group).
- **`FOURIER`**: (Shim module) For `GetFFTComplexDims` and `GetFFTDims`. These routines from the shim module internally call the `_NEW` versions from `FOURIER_NEW`.
- **`FOURIER_NEW`**: This is the actual FFT interface module. `InitializeFFT` calls `PlanFFT_NEW` from this module to set up FFTW3 plans, using `FFT_STD_STATE` (a configuration object from `FOURIER_NEW`) and planning flags like `FFTW3_ONEPLAN_ESTIMATE`.
- **`OutputFiles`**: For `outputUnit` (Fortran I/O unit for standard output).
- **`CellInfo`**: For cell lattice vector lengths (`cell%lengthA`, `cell%lengthB`, `cell%lengthC`) used in `SizeSystem`.

The `InitializeFFT` routine is a critical setup step called early in the program (e.g., from `Initializer` module's `SetupAllModules`) to prepare for all subsequent FFT operations. The choice of `dimType` and `energyCutoff`/`gridSpacing` will determine the size and efficiency of FFTs throughout the simulation.
