# Overview

The `Fourier_NEW` module (in file `Fourier_v1.f90`) provides a Fortran interface to the FFTW (Fastest Fourier Transform in the World) version 3.x library. This specific version of the module is designed for **serial execution**. It handles the planning of Fourier transforms, their execution, and the subsequent cleanup of resources. The module aims to provide a structured way to perform FFTs for various quantities within the application.

The general workflow involves:
1.  Planning the FFT for specific dimensions and types using `PlanFFT_NEW`.
2.  Executing the planned transform as needed using the generic `FFT_NEW` interface, which dispatches to the correct underlying 3D or 4D, forward or backward transform subroutines.
3.  Cleaning up (releasing) the FFT plan resources using `CleanFFT_NEW` when they are no longer needed.

# Key Components

- **`MODULE Fourier_NEW`**: The main module encapsulating all FFT functionalities for serial operations.

- **`TYPE FFT_CONFIG`**:
  - **Description:** A derived data type that stores the configuration and state for a particular FFT setup. This includes dimensions, planning mode, normalization constants, and FFTW plan handles.
  - **Fields:**
    - `fftAlgo :: INTEGER`: The algorithm/mode flag used for creating FFTW plans (e.g., `FFTW_ESTIMATE`, `FFTW_MEASURE`).
    - `totalDimX, totalDimY, totalDimZ :: INTEGER`: The global dimensions of the 3D real-space grid.
    - `localDimZ, localDimZOff :: INTEGER`: For this serial version, `localDimZ` will be equal to `totalDimZ`, and `localDimZOff` will be 0.
    - `iCountFFT :: INTEGER`: A counter, possibly for tracking the number of FFTs performed.
    - `offset :: INTEGER`: Stores the parity (`MOD(totalDimX, 2)`) of the X-dimension, used to correctly determine array sizes for real-to-complex transforms.
    - `normConst`: The normalization constant applied to the FFT results. Its type (`INTEGER(KIND=DP)` or `REAL(KIND=DP)`) depends on the `FFT_NORM_EXACT` preprocessor macro.
    - `fftw3GlobFWD, fftw3GlobBWD :: INTEGER(kind=8)`: Handles (pointers stored as integers) to the pre-computed FFTW plans for forward (real-to-complex) and backward (complex-to-real) transforms.

- **`FFT_STD_STATE :: TYPE(FFT_CONFIG), SAVE`**:
  - **Description:** A module-level variable of `FFT_CONFIG` type, saved across calls. It likely holds the configuration for a standard, commonly used FFT plan within the application.

- **`INTERFACE FFT_NEW`**:
  - **Description:** A generic public interface for performing FFTs. Client code calls `FFT_NEW(config, array, transform)`, and the interface selects the appropriate underlying subroutine based on the rank and type of `array` and `transform`.
  - **Module Procedures:** (These are implemented as subroutines)
    - `ForwardFFT_4D`
    - `BackFFT_4D`
    - `ForwardFFT_3D`
    - `BackFFT_3D`

- **`SUBROUTINE PlanFFT_NEW(config, mode, dimX, dimY, dimZ)`**:
  - **Description:** Initializes an `FFT_CONFIG` instance by creating FFTW plans for the given dimensions and planning mode. For "one-plan" modes, it allocates temporary arrays to create the plans.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`: The configuration object to be initialized.
    - `mode :: INTEGER, INTENT(IN)`: The FFTW planning mode (e.g., `FFTW3_ONEPLAN_ESTIMATE`).
    - `dimX, dimY, dimZ :: INTEGER, INTENT(IN)`: Global dimensions of the real-space grid.

- **`SUBROUTINE GetFFTDims_NEW(config, dimX, dimY, dimZ, zoff)`**:
  - **Description:** Retrieves the real-space grid dimensions from an `FFT_CONFIG` object. For serial, `dimZ` is global and `zoff` is 0.
  - **Arguments (Output):**
    - `dimX, dimY, dimZ, zoff :: INTEGER`

- **`SUBROUTINE GetFFTComplexDims_NEW(config, dimX, dimY, locDimZ, locDimZOffset)`**:
  - **Description:** Retrieves reciprocal-space grid dimensions. `dimX` is `totalDimX/2 + 1`. For serial, `locDimZ` is global and `locDimZOffset` is 0.
  - **Arguments (Output):**
    - `dimX, dimY, locDimZ, locDimZOffset :: INTEGER`

- **`SUBROUTINE ForwardFFT_4D(config, array, transform)`** (Private):
  - **Description:** Implements 4D forward FFT by calling `ForwardFFT_3D` for each slice along the 4th dimension (typically spin).
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN)`
    - `transform :: COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), SIZE(array,3), SIZE(array,4)), INTENT(OUT)`

- **`SUBROUTINE BackFFT_4D(config, array, transform)`** (Private):
  - **Description:** Implements 4D backward FFT by calling `BackFFT_3D` for each slice along the 4th dimension.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: COMPLEX(kind=DP), DIMENSION(:,:,:,:), INTENT(IN)`
    - `transform :: REAL(kind=DP), DIMENSION(2*(SIZE(array,1)-1)+config%offset, SIZE(array,2), SIZE(array,3), SIZE(array,4)), INTENT(OUT)`

- **`SUBROUTINE ForwardFFT_3D(config, array, transform)`** (Private):
  - **Description:** Executes a 3D forward (real-to-complex) FFT. If a global plan exists in `config` (from "one-plan" modes), it's used. Otherwise, a temporary plan is created, executed, and destroyed. Applies normalization.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)`
    - `transform :: COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), SIZE(array,3)), INTENT(OUT)`

- **`SUBROUTINE BackFFT_3D(config, array, transform)`** (Private):
  - **Description:** Executes a 3D backward (complex-to-real) FFT. Similar to `ForwardFFT_3D`, it uses either a global plan from `config` or a temporary plan.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(INOUT)`
    - `array :: COMPLEX(kind=DP), DIMENSION(:,:,:), INTENT(IN)`
    - `transform :: REAL(kind=DP), DIMENSION(2*(SIZE(array,1)-1)+config%offset, SIZE(array,2), SIZE(array,3)), INTENT(OUT)`

- **`SUBROUTINE CleanFFT_NEW(config)`**:
  - **Description:** Destroys the global FFTW plans stored in `config` if they were created using "one-plan" modes.
  - **Arguments:**
    - `config :: TYPE(FFT_CONFIG), INTENT(IN)`

- **`SUBROUTINE ImportWisdom(ioNum)` / `SUBROUTINE ExportWisdom(ioNum)`**:
  - **Description:** Placeholder subroutines for FFTW wisdom import/export functionality (currently marked as TODO).

# Important Variables/Constants

- **`FFTW3_ESTIMATE, FFTW3_ONEPLAN_ESTIMATE, FFTW3_ONEPLAN_MEASURE, FFTW3_ONEPLAN_PATIENT :: INTEGER, PARAMETER`**: Constants defining different planning strategies for FFTW. `FFTW_ESTIMATE` is faster to plan but may result in slower execution. `FFTW_MEASURE` and `FFTW_PATIENT` take longer to plan but can yield faster transforms. The "ONEPLAN" versions create persistent plans stored in `FFT_CONFIG`.
- **`FFT_NORM_EXACT`**: A preprocessor macro. If defined, `normConst` is an integer and normalization is by division. Otherwise, `normConst` is real and normalization is by multiplication (inverse factor).
- **`DEBUG_FFT_INTER`**: A preprocessor macro. If defined, enables printing of debug information about FFT plans.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- **FFTW3 Library**: This module is a direct wrapper for the serial FFTW3 library.
    - `fftw3.f`: The Fortran include file for FFTW3.
    - It uses various `dfftw_*` routines from the library, such as `dfftw_plan_dft_r2c_3d`, `dfftw_plan_dft_c2r_3d`, `dfftw_execute_dft_r2c`, `dfftw_execute_dft_c2r`, and `dfftw_destroy_plan`.
- **Timer/Output**: Calls `StartClock` and `StopClock` (likely from a `Timer` or `Output` module) for profiling FFT execution times.

This module provides the core FFT capabilities for the application when running in serial mode. It contrasts with `Fourier_para_v1.f90` which handles parallel FFTs using FFTW3's MPI interface.
