# Overview

The `MPI_Functions` module centralizes MPI (Message Passing Interface) operations required for parallel execution of the PROFESS code. It provides a set of wrapper subroutines for common MPI tasks, such as initialization, finalization, data broadcasting, and global reductions.

A key feature is the setup of two MPI communicators:
1.  `MPI_COMM_WORLD`: The default communicator encompassing all processors.
2.  `MPI_GFFT_WORLD`: A potentially smaller communicator, derived by splitting `MPI_COMM_WORLD`. This group is intended for processors that will participate in global Fast Fourier Transform (FFT) calculations, allowing for a flexible allocation of resources if not all processors are needed for this specific task.

The actual MPI calls are conditionally compiled using the `__USE_PARALLEL` preprocessor macro. If this macro is not defined, the code generally runs in serial mode, and these MPI wrapper routines will have no effect or issue a warning.

# Key Components

- **`MODULE MPI_Functions`**: The main container for MPI-related utilities.

- **`SUBROUTINE InitializeMPI()`**:
  - **Description:** Initializes the MPI environment. It calls `MPI_INIT`, determines the rank (`rankGlobal`) and size (`sizeGlobal`) of the `MPI_COMM_WORLD`. If `numProcGFFT` (number of processors for global FFT group) is specified and positive, it splits `MPI_COMM_WORLD` to create `MPI_GFFT_WORLD` and determines the rank (`rankGFFT`) and size (`sizeGFFT`) within this new communicator. If `numProcGFFT` is -1, `MPI_GFFT_WORLD` becomes equivalent to `MPI_COMM_WORLD`. It also sets the `outRank` variable (from `OutputFiles` module) to `rankGlobal`.
  - **Side Effects:** Initializes MPI, sets `rankGlobal`, `sizeGlobal`, `MPI_GFFT_WORLD`, `rankGFFT`, `sizeGFFT`, `indGFFT`, and `outRank`.

- **`SUBROUTINE QuitMPI()`**:
  - **Description:** Finalizes the MPI environment by calling `MPI_FINALIZE`. This should be called before the program exits.

- **`SUBROUTINE BcastReal(value, valueLen)`**:
  - **Description:** Broadcasts a 1D array of real numbers (`REAL(KIND=DP)`) from the root processor (rank 0 in `MPI_COMM_WORLD`) to all other processors in `MPI_COMM_WORLD`.
  - **Arguments:**
    - `value :: REAL(KIND=DP), DIMENSION(valueLen), INTENT(INOUT)`: The array to be broadcast/received.
    - `valueLen :: INTEGER, INTENT(IN)`: The number of elements in the array.

- **`SUBROUTINE BcastReal_Dim2(value, valueLen)`**:
  - **Description:** Broadcasts a 2D real array, specifically dimensioned as `(3, valueLen)`, from the root to all other processors.
  - **Arguments:**
    - `value :: REAL(KIND=DP), DIMENSION(3,valueLen), INTENT(INOUT)`.
    - `valueLen :: INTEGER, INTENT(IN)`: The size of the second dimension.

- **`SUBROUTINE BcastCharacter(word, wordLen)`**:
  - **Description:** Broadcasts a character array (string) from the root to all other processors.
  - **Arguments:**
    - `word :: CHARACTER, DIMENSION(wordLen), INTENT(INOUT)`: The character array.
    - `wordLen :: INTEGER, INTENT(IN)`: The length of the character array.

- **`INTERFACE ReduceRealLevel1`**:
  - **Description:** A generic interface for performing sum reductions on real-valued data across all processors in `MPI_COMM_WORLD`. The result is available on all processors (due to `MPI_ALLREDUCE`).
  - **Module Procedures:**
    - `ReduceReal_single(value)`: For a single scalar `REAL(KIND=DP)`.
    - `ReduceReal_array(value, nvalue)`: For a 1D `REAL(KIND=DP)` array.
    - `ReduceRealLevel1_2Darray(value, nx, ny)`: For a 2D `REAL(KIND=DP)` array.
    - `ReduceRealLevel1_3Darray(value, nx, ny, nz)`: For a 3D `REAL(KIND=DP)` array.

- **`SUBROUTINE ReduceInteger_single(value)`**:
  - **Description:** Performs a sum reduction on a single integer value across all processors in `MPI_COMM_WORLD`. The result is available on all processors.
  - **Arguments:** `value :: INTEGER, INTENT(INOUT)`.

# Important Variables/Constants

- **Module-Level MPI Information:**
    - `rankGlobal :: INTEGER`: Rank of the current MPI process within `MPI_COMM_WORLD` (0 to `sizeGlobal-1`). Initialized to 0.
    - `sizeGlobal :: INTEGER`: Total number of MPI processes in `MPI_COMM_WORLD`. Initialized to 1.
    - `MPI_GFFT_WORLD :: INTEGER`: MPI Communicator handle for the (potentially smaller) group of processors designated for global FFT operations.
    - `rankGFFT :: INTEGER`: Rank of the current MPI process within `MPI_GFFT_WORLD`. Initialized to 0.
    - `sizeGFFT :: INTEGER`: Total number of MPI processes in `MPI_GFFT_WORLD`. Initialized to 1.
    - `indGFFT :: INTEGER`: The color or group index used when splitting `MPI_COMM_WORLD` to create `MPI_GFFT_WORLD`.
    - `numProcGFFT :: INTEGER`: User-defined number of processors to be included in the `MPI_GFFT_WORLD` communicator. A value of -1 (default) implies `MPI_GFFT_WORLD` will be the same as `MPI_COMM_WORLD`.
    - `mpiErr :: INTEGER`: Stores the error status returned by MPI calls.
- **Control Variables:**
    - `warning :: INTEGER`: A flag (default 0). If greater than 0 and `DEBUG_PROFESS` is defined, MPI functions will print a warning if called in a serial build (where `__USE_PARALLEL` is not defined).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- **`OutputFiles`**: For `errorUnit` (Fortran I/O unit for error messages) and for setting the `outRank` variable (which is declared in `OutputFiles` but set by `InitializeMPI`).
- **MPI Library**: This module is fundamentally a wrapper around an MPI library.
    - It includes `mpif.h` when `__USE_PARALLEL` is defined.
    - It uses MPI routines such as `MPI_INIT`, `MPI_FINALIZE`, `MPI_COMM_RANK`, `MPI_COMM_SIZE`, `MPI_BCAST`, `MPI_ALLREDUCE`, and `MPI_COMM_SPLIT`.
- **Preprocessor Macros**:
    - `__USE_PARALLEL`: If defined, the MPI calls are compiled and executed. Otherwise, the MPI-related code is skipped.
    - `DEBUG_PROFESS`: If defined (and `__USE_PARALLEL` is not), enables warnings when MPI wrapper functions are called in a serial build, provided `warning > 0`.

This module is essential for enabling parallel execution of the PROFESS code. Other modules that need to perform parallel communication or synchronization will `USE MPI_Functions` and call these wrapper routines. The main program will typically call `InitializeMPI` at the very beginning and `QuitMPI` at the very end.
