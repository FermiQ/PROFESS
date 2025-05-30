# Overview

The `MathFunctions` module provides a collection of common mathematical operations and utility functions frequently used in scientific computing, particularly within the context of materials simulations. These include vector and matrix algebra for 3D systems (cross product, norm, determinant, inverse for both real and complex types), specialized physics-related functions like the Lindhard response function (`LindG`) and its derivative (`GPrime`), and general utilities such as string manipulation (`Uppercase`), array operations (`MinMaxVal`, `stepfun`), metric tensor calculation (`Metric3d`), matrix diagonalization (`Diag`), and matrix symmetrization (`Symmetrize`).

Many functions are provided via generic interfaces (e.g., `Cross`, `Norm`, `Det`, `Inverse`) that dispatch to the appropriate real or complex version based on argument types.

# Key Components

- **`MODULE MathFunctions`**: The main container for all mathematical functions and subroutines.

- **Generic Interfaces for Vector/Matrix Operations:**
    - **`INTERFACE Cross`**: Computes the 3D vector cross product.
        - `FUNCTION CrossProductComplex(a,b)`: For complex 3-vectors.
        - `FUNCTION CrossProductReal(a,b)`: For real 3-vectors.
    - **`INTERFACE Norm`**: Computes the Euclidean norm (length) of a 3D vector.
        - `FUNCTION NormComplex(a)`: For complex 3-vectors (`sqrt(a .dot. a)`).
        - `FUNCTION NormReal(a)`: For real 3-vectors (`sqrt(a .dot. a)`).
    - **`INTERFACE Det`**: Computes the determinant of a 3x3 matrix.
        - `FUNCTION DetComplex(M)`: For complex 3x3 matrices.
        - `FUNCTION DetReal(M)`: For real 3x3 matrices.
    - **`INTERFACE Inverse`**: Computes the inverse of a 3x3 matrix.
        - `FUNCTION InverseComplex(M)`: For complex 3x3 matrices.
        - `FUNCTION InverseReal(M)`: For real 3x3 matrices.

- **Specialized Physics/Math Functions:**
    - **`FUNCTION LindG(eta, lambda, mu) RESULT(value)`**:
      - **Description:** Calculates the Lindhard response function `G(eta)` or a similar kernel form, often used in Wang-Teter KEDF. It includes specific treatments for `eta` near 0, near 1 (singularity), and a Taylor expansion for large `eta` (>3.65). The parameters `lambda` and `mu` are typically related to TF and vW coefficients.
      - **Arguments:** `eta`, `lambda`, `mu` (all `REAL(KIND=DP)`).
    - **`FUNCTION GPrime(eta, mu) RESULT(value)`**:
      - **Description:** Calculates the first derivative of the `LindG` function with respect to `eta`.
      - **Arguments:** `eta`, `mu` (both `REAL(KIND=DP)`).
    - **`FUNCTION Volume(M) RESULT(value)`**:
      - **Description:** Calculates the volume of a parallelepiped defined by the column (or row) vectors of a 3x3 real matrix `M` (typically the cell lattice vectors). Calculated as `abs(Det(M))`.
      - **Arguments:** `M :: REAL(KIND=DP), DIMENSION(3,3), INTENT(IN)`.

- **Utility Subroutines and Functions:**
    - **`SUBROUTINE Uppercase(text)`**: Converts an input character string `text` to all uppercase letters in place.
    - **`FUNCTION Vecmul(matrix, vector) RESULT(vector_out)`**: Multiplies a 3x3 real matrix by a 3-component real vector.
    - **`FUNCTION MinMaxVal(data, flag, useData) RESULT(value)`**:
      - **Description:** Finds the global minimum, maximum, or sum of elements in a 3D real array `data`. If `__USE_PARALLEL` is defined, it performs an MPI reduction (`MPI_MIN`, `MPI_MAX`, or `MPI_SUM`) to get the global value across all processors. An optional logical array `useData` can mask elements from the calculation.
      - **Arguments:** `data(:,:,:)`, `flag` (CHARACTER: 'MIN', 'MAX', or 'SUM'), `useData(:,:,:)` (OPTIONAL).
    - **`SUBROUTINE Metric3d(cellReal, rmet, gmet)`**:
      - **Description:** Computes the real-space metric tensor `rmet_ij = a_i .dot. a_j` and the reciprocal-space metric tensor `gmet_ij = b_i .dot. b_j`, where `a_i` are real-space lattice vectors (columns of `cellReal`) and `b_i` are reciprocal lattice vectors.
      - **Arguments:** `cellReal(3,3)` (IN), `rmet(3,3)` (OUT), `gmet(3,3)` (OUT).
    - **`SUBROUTINE Diag(n, a, eig, eigv)`**:
      - **Description:** Computes all eigenvalues (`eig`) and eigenvectors (`eigv`, stored column-wise) of a real symmetric `n x n` matrix `a`. It uses the LAPACK routine `DSPEV`. The input matrix `a` is provided in packed upper triangular form to `DSPEV`.
      - **Arguments:** `n` (IN), `a(n,n)` (IN), `eig(n)` (OUT), `eigv(n,n)` (OUT).
    - **`SUBROUTINE Symmetrize(n, a)`**:
      - **Description:** Symmetrizes an `n x n` real matrix `a` in place by setting `a_ij = a_ji = (old_a_ij + old_a_ji)/2`.
      - **Arguments:** `n` (IN), `a(n,n)` (INOUT).
    - **`FUNCTION stepfun(array) RESULT(step_array)`**:
      - **Description:** An element-wise step function. For each element in the input 3D real `array`, the corresponding element in `step_array` is 1.0 if `array(i,j,k) >= 0`, and 0.0 otherwise.
      - **Arguments:** `array(:,:,:) :: REAL(KIND=DP), INTENT(IN)`.

# Important Variables/Constants

This module primarily defines functions and subroutines. It imports `DP` and `pi` from the `CONSTANTS` module.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter) and `pi`.
- **`TIMER`**: For `stopWatch` type, `TimerStart`, and `TimerStop` (used for profiling, though not shown in all functions).
- **LAPACK (External Library)**: The `Diag` subroutine relies on the LAPACK routine `DSPEV` for matrix diagonalization. This implies a need to link against a LAPACK library.
- **MPI (External Library)**: The `MinMaxVal` function uses MPI routines (`MPI_ALLREDUCE`, etc.) for global reductions when the `__USE_PARALLEL` preprocessor macro is defined. This requires linking against an MPI library. It includes `mpif.h` if `__USE_PARALLEL` is defined.

This module provides a suite of mathematical tools used throughout the PROFESS codebase for various calculations involving vectors, matrices, and specific mathematical forms relevant to DFT and KEDF theory.
