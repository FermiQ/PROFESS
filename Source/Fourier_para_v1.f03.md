# Overview

The `Fourier_new` module provides a Fortran interface to the FFTW (Fastest Fourier Transform in the World) version 3.x library, specifically tailored for parallel execution using MPI. It encapsulates the necessary steps for performing FFTs: planning the transform, executing it, and cleaning up resources. This module is intended to replace older FFT implementations and offers a more modern, state-free approach (though it uses a saved `FFT_CONFIG` type, `FFT_STD_STATE`, for a default plan).

The module defines a configuration type `FFT_CONFIG` to hold plan details, dimensions, and workspace pointers. It supports forward (real to complex) and backward (complex to real) transforms for 3D and 4D data (where the 4th dimension is typically spin).

# Key Components

- **`MODULE Fourier_new`**: The main module.

- **`TYPE FFT_CONFIG`**:
  - **Description:** A derived data type to store all information related to a specific FFT setup. This includes dimensions, MPI configuration, FFTW plans, and pointers to data buffers.
  - **Fields:**
    - `fftAlgo :: INTEGER`: Specifies the FFTW planning mode (e.g., `FFTW_ESTIMATE`, `FFTW_MEASURE`).
    - `totalDimX, totalDimY, totalDimZ :: INTEGER(C_INTPTR_T)`: Global dimensions of the 3D real grid.
    - `localDimZ, localDimZOff :: INTEGER(C_INTPTR_T)`: The size of the local slab in the Z-dimension and its offset for the current MPI process.
    - `iCountFFT :: INTEGER`: A counter for FFT operations, likely for debugging or statistics.
    - `normConst`: Normalization constant for the transforms. Can be `INTEGER(KIND=DP)` or `REAL(KIND=DP)` based on `FFT_NORM_EXACT` macro.
    - `fftw3GlobFWD, fftw3GlobBWD :: TYPE(C_PTR)`: C pointers to the FFTW plans for forward (real-to-complex) and backward (complex-to-real) transforms, respectively.
    - `cdata :: TYPE(C_PTR)`: C pointer to the contiguous memory block allocated by FFTW for the transform data (used for in-place transforms).
    - `in :: REAL(C_DOUBLE), POINTER`: Fortran pointer aliased to `cdata`, interpreted as a 3D real array.
    - `out :: COMPLEX(C_DOUBLE_COMPLEX), POINTER`: Fortran pointer aliased to `cdata`, interpreted as a 3D complex array.

- **`FFT_STD_STATE :: TYPE(FFT_CONFIG), SAVE`**:
  - **Description:** A module-level variable of type `FFT_CONFIG` with the `SAVE` attribute. This likely holds the plan and configuration for a standard, frequently used FFT setup throughout the application.

- **`INTERFACE FFT_NEW`**:
  - **Description:** A generic public interface for performing FFTs. It allows users to call `FFT_NEW(config, array, transform)` where `config` is an `FFT_CONFIG` instance, and the specific transform (3D/4D, forward/backward) is chosen based on the types and ranks of `array` and `transform`.
  - **Module Procedures:** (These are now subroutines, not functions as in the older `Fourier` module)
    - `ForwardFFT_4D`
    - `BackFFT_4D`
    - `ForwardFFT_3D`
    - `BackFFT_3D`

- **`SUBROUTINE PlanFFT_NEW(MPI_COMM, config, mode, dimX, dimY, dimZ)`**:
  - **Description:** Initializes an `FFT_CONFIG` instance (`config`) by creating FFTW plans for a given set of global dimensions (`dimX, dimY, dimZ`) and MPI communicator (`MPI_COMM`). It determines local data distribution and allocates memory.
  - **Arguments:**
    - `MPI_COMM :: INTEGER, INTENT(IN)`: The MPI communicator for parallel transforms.
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`: The FFT configuration to be initialized.
    - `mode :: INTEGER, INTENT(IN)`: FFTW planning strategy (e.g., `FFTW3_ONEPLAN_MEASURE`).
    - `dimX, dimY, dimZ :: INTEGER, INTENT(IN)`: Global dimensions of the real-space grid.

- **`SUBROUTINE GetFFTDims_NEW(config, dimX, dimY, dimZ, zoff)`**:
  - **Description:** Retrieves the real-space grid dimensions (global `dimX`, `dimY` and local `dimZ`, `zoff`) from a given `FFT_CONFIG` instance.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(IN)`
    - `dimX, dimY, dimZ, zoff :: INTEGER, INTENT(OUT)`

- **`SUBROUTINE GetFFTComplexDims_NEW(config, dimX, dimY, locDimZ, locDimZOffset)`**:
  - **Description:** Retrieves the reciprocal-space grid dimensions from an `FFT_CONFIG` instance. `dimX` will be `totalDimX/2 + 1`.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(IN)`
    - `dimX, dimY, locDimZ, locDimZOffset :: INTEGER, INTENT(OUT)`

- **`SUBROUTINE ForwardFFT_4D(config, array, transform)`** (Private):
  - **Description:** Performs a 4D forward FFT by iterating over the 4th dimension (spin) and calling `ForwardFFT_3D` for each 3D slice.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN)`
    - `transform :: COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), SIZE(array,3), SIZE(array,4)), INTENT(OUT)`

- **`SUBROUTINE BackFFT_4D(config, array, transform)`** (Private):
  - **Description:** Performs a 4D backward FFT by iterating over the 4th dimension and calling `BackFFT_3D`.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: COMPLEX(kind=DP), DIMENSION(:,:,:,:), INTENT(IN)`
    - `transform :: REAL(kind=DP), DIMENSION(config%totalDimX,config%totalDimY,config%localDimZ, SIZE(array,4)), INTENT(OUT)`

- **`SUBROUTINE ForwardFFT_3D(config, array, transform)`** (Private):
  - **Description:** Executes a 3D forward (real-to-complex) FFT using the pre-established plan in `config`. Data is copied into the FFTW buffer, transformed in-place, and then copied to the output `transform` array with normalization.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)`
    - `transform :: COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), SIZE(array,3)), INTENT(OUT)`

- **`SUBROUTINE BackFFT_3D(config, array, transform)`** (Private):
  - **Description:** Executes a 3D backward (complex-to-real) FFT. Data is copied to the FFTW buffer, transformed in-place, and then copied to `transform`. Normalization for backward transform is typically handled by FFTW itself or implicitly by the forward transform's normalization.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: COMPLEX(kind=DP), DIMENSION(:,:,:), INTENT(IN)`
    - `transform :: REAL(kind=DP), DIMENSION(config%totalDimX,config%totalDimY, config%localDimZ), INTENT(OUT)`

- **`SUBROUTINE CleanFFT_NEW(config)`**:
  - **Description:** Releases resources associated with an `FFT_CONFIG` instance, including destroying FFTW plans and freeing allocated memory.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(IN)`

- **`SUBROUTINE ImportWisdom(ioNum)` / `SUBROUTINE ExportWisdom(ioNum)`**:
  - **Description:** Stubs for importing/exporting FFTW "wisdom" (saved information about optimal plans), marked as TODO.

# Important Variables/Constants

- **`FFTW3_ONEPLAN_ESTIMATE, FFTW3_ONEPLAN_MEASURE, FFTW3_ONEPLAN_PATIENT :: INTEGER, PARAMETER`**: Constants representing different FFTW planning modes, trading off planning time for runtime efficiency.
- **`FFT_NORM_EXACT`**: A preprocessor macro that, if defined, changes how the normalization constant `normConst` in `FFT_CONFIG` is calculated and applied (integer vs. real division).
- **`DEBUG_FFT_INTER`**: A preprocessor macro that, if defined, enables debug print statements.
- **`PROFESS`**: A preprocessor macro that, if defined, enables calls to `StartClock` and `StopClock` (profiling).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- **`iso_c_binding`**: Intrinsic Fortran module for interoperability with C, necessary for interfacing with FFTW.
- **FFTW3 Library**: This module is fundamentally a wrapper around the FFTW3 library, specifically its MPI-enabled version. It uses:
    - `fftw3-mpi.f03`: Fortran include file for FFTW3 MPI.
    - Various `fftw_mpi_*` and `fftw_*` functions (e.g., `fftw_mpi_local_size_3d`, `fftw_alloc_complex`, `fftw_mpi_plan_dft_r2c_3d`, `fftw_mpi_execute_dft_r2c`, `fftw_destroy_plan`, `fftw_free`).
- **MPI**: Requires an MPI library for parallel execution.
    - `mpif.h`: MPI Fortran include file (used if `__USE_PARALLEL` is defined, though this macro isn't explicitly used in the provided code for conditional MPI calls, `fftw3-mpi.f03` implies MPI usage).
- **Timer/Output (Conditional)**: If `PROFESS` macro is defined, it calls `StartClock` and `StopClock`, which are assumed to be from a timing/profiling utility module (e.g., `Timer` or `Output`).

The typical workflow with this module for a given `FFT_CONFIG` instance (e.g., `FFT_STD_STATE`):
1. Call `PlanFFT_NEW` once at the beginning to create and store FFTW plans in the `FFT_CONFIG` instance.
2. Use `GetFFTDims_NEW` and `GetFFTComplexDims_NEW` if grid dimension information is needed.
3. Call `FFT_NEW(config, input_array, output_array)` to perform transforms. The interface dispatches to the appropriate 3D/4D, forward/backward subroutine.
4. Call `CleanFFT_NEW` at the end of the program to release FFTW resources associated with that `FFT_CONFIG` instance.
