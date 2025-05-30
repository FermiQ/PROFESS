# Overview

The `Fourier` module serves as a compatibility layer or "shim" for performing Fast Fourier Transforms (FFTs). Its primary purpose is to provide a consistent interface for FFT operations while internally utilizing a newer, state-free FFT implementation found in the `FOURIER_NEW` module. The comment "This is an ugly shim" suggests it's intended as a temporary solution before client code is updated to use `FOURIER_NEW` directly.

The module handles planning of FFTs, retrieving FFT grid dimensions, and performing forward (real to complex) and backward (complex to real) transforms for both 3D and 4D arrays (where the 4th dimension is typically spin).

# Key Components

- **`MODULE Fourier`**: The main container module.

- **`INTERFACE FFT`**:
  - **Description:** Provides a generic public interface for FFT operations. Based on the type of the input array (real or complex), it automatically selects the appropriate forward or backward transform. This simplifies FFT calls in the client code.
  - **Module Procedures:**
    - `ForwardFFT_4D`: For transforming 4D real arrays to 4D complex arrays.
    - `BackFFT_4D`: For transforming 4D complex arrays back to 4D real arrays.
    - `ForwardFFT_3D`: For transforming 3D real arrays to 3D complex arrays.
    - `BackFFT_3D`: For transforming 3D complex arrays back to 3D real arrays.

- **`SUBROUTINE PlanFFT(MPI_COMM, dimX, dimY, dimZ)`** (Parallel version) / **`SUBROUTINE PlanFFT(dimX, dimY, dimZ)`** (Serial version):
  - **Description:** Initializes the FFT plan(s) required by the underlying FFT library (likely FFTW via `FOURIER_NEW`). This typically involves analyzing the transform size and type to choose an optimal algorithm. It calls `PlanFFT_NEW` from the `FOURIER_NEW` module.
  - **Arguments:**
    - `MPI_COMM :: INTEGER, INTENT(IN)`: MPI communicator (for parallel version).
    - `dimX, dimY, dimZ :: INTEGER, INTENT(IN)`: The global dimensions of the real-space grid.

- **`SUBROUTINE GetFFTDims(dimX, dimY, locDimZ, locDimZOffset)`**:
  - **Description:** Retrieves the dimensions of the real-space FFT grid, including local dimensions (`locDimZ`) and offset (`locDimZOffset`) for the current MPI process. It calls `GetFFTDims_NEW`.
  - **Arguments (Output):**
    - `dimX, dimY :: INTEGER`: Global dimensions of the real-space grid.
    - `locDimZ, locDimZOffset :: INTEGER`: Local z-dimension and its offset for the current process.

- **`SUBROUTINE GetFFTComplexDims(dimX, dimY, locDimZ, locDimZOffset)`**:
  - **Description:** Retrieves the dimensions of the reciprocal-space (complex) FFT grid, including local dimensions and offset. It calls `GetFFTComplexDims_NEW`.
  - **Arguments (Output):**
    - `dimX, dimY :: INTEGER`: Dimensions of the complex grid (note: `dimX` is typically `global_dimX/2 + 1`).
    - `locDimZ, locDimZOffset :: INTEGER`: Local z-dimension and its offset for the complex grid.

- **`FUNCTION ForwardFFT_4D(array) RESULT(transform)`** (Private):
  - **Description:** Performs a forward FFT on a 4D real array (e.g., real-space density with spin channels) to produce a 4D complex array. It iterates over the 4th dimension (spin) and calls `FFT_NEW` for each 3D slice. Marked for deprecation.
  - **Arguments:**
    - `array :: REAL(kind=DP), DIMENSION(:,:,:,:)`: Input 4D real array.
  - **Return Value:**
    - `transform :: COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), SIZE(array,3), SIZE(array,4))`: Output 4D complex array.

- **`FUNCTION BackFFT_4D(array) RESULT(transform)`** (Private):
  - **Description:** Performs a backward FFT on a 4D complex array to produce a 4D real array. Iterates over the 4th dimension and calls `FFT_NEW`. Marked for deprecation.
  - **Arguments:**
    - `array :: COMPLEX(kind=DP), DIMENSION(:,:,:,:)`: Input 4D complex array.
  - **Return Value:**
    - `transform :: REAL(kind=DP), DIMENSION(FFT_STD_STATE%totalDimX, FFT_STD_STATE%totalDimY, FFT_STD_STATE%localDimZ, SIZE(array,4))`: Output 4D real array.

- **`FUNCTION ForwardFFT_3D(array) RESULT(transform)`** (Private):
  - **Description:** Performs a forward FFT on a 3D real array. Calls `FFT_NEW`.
  - **Arguments:**
    - `array :: REAL(kind=DP), DIMENSION(:,:,:)`: Input 3D real array.
  - **Return Value:**
    - `transform :: COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), SIZE(array,3))`: Output 3D complex array.

- **`FUNCTION BackFFT_3D(array) RESULT(transform)`** (Private):
  - **Description:** Performs a backward FFT on a 3D complex array. Calls `FFT_NEW`.
  - **Arguments:**
    - `array :: COMPLEX(kind=DP), DIMENSION(:,:,:)`: Input 3D complex array.
  - **Return Value:**
    - `transform :: REAL(kind=DP), DIMENSION(FFT_STD_STATE%totalDimX, FFT_STD_STATE%totalDimY, FFT_STD_STATE%localDimZ)`: Output 3D real array.

- **`SUBROUTINE CleanFFT`**:
  - **Description:** Frees any memory or resources associated with the FFT plans. Calls `CleanFFT_NEW`.

# Important Variables/Constants

- **`offset :: INTEGER, SAVE`**: Stores `MOD(dimX, 2)` after `PlanFFT` is called. This is likely used to handle the ambiguity in real-space grid size (2k or 2k+1) when the reciprocal-space size is k+1.
- **`iCountFFT :: INTEGER, SAVE`**: A counter that is incremented each time `ForwardFFT_3D` or `BackFFT_3D` is called. Its exact purpose is not detailed but might be for debugging or tracking FFT usage.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- **`FOURIER_NEW`**: This is the primary dependency. The `Fourier` module is essentially a wrapper around routines from `FOURIER_NEW`. Specific routines used include:
    - `PlanFFT_NEW`
    - `GetFFTDims_NEW`
    - `GetFFTComplexDims_NEW`
    - `FFT_NEW` (the core transform routine)
    - `CleanFFT_NEW`
  It also uses named constants/types from `FOURIER_NEW` like `FFT_STD_STATE` (likely a derived type holding FFT plan state) and `FFTW3_ONEPLAN_MEASURE` (likely a flag for FFTW planning).

- **Parallelism**: The `PlanFFT` subroutine has distinct versions for serial and parallel execution, controlled by the `__USE_PARALLEL` preprocessor macro. The parallel version takes an `MPI_COMM` argument.

The typical workflow with this module would be:
1. Call `PlanFFT` once at the beginning of a calculation to initialize FFT plans.
2. Use `GetFFTDims` and `GetFFTComplexDims` if needed to ascertain array dimensions.
3. Call `FFT(array)` to perform transforms. The interface will dispatch to the correct underlying 3D/4D, forward/backward function.
4. Call `CleanFFT` at the very end of the program to release resources.
