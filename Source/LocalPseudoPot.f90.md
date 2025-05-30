# Overview

The `LocalPseudoPot` module provides functionalities for evaluating numerically defined local pseudopotentials, primarily in reciprocal space. A local pseudopotential `V_local(r)` (or `V_local(q)` in reciprocal space) describes the interaction of an electron with an ion core. This module allows the lookup of pseudopotential values `V(q)` and their first derivatives `dV(q)/dq` at arbitrary `q` magnitudes using interpolation.

A key aspect of the reciprocal space lookup is the handling of the Coulombic tail (`-4*pi*Z/q^2`). To enable accurate spline interpolation near `q=0` where `V(q)` diverges, the module works with a modified quantity `V(q) + 4*pi*Z/q^2`, which is finite at `q=0`. After interpolation of this modified quantity, the Coulomb tail is subtracted back to yield the actual `V(q)`.

# Key Components

- **`MODULE LocalPseudoPot`**: The main container for pseudopotential lookup functions.

- **`FUNCTION PseudoPotLookup(psp, lookUp) RESULT(value)`**:
  - **Description:** This is the primary public interface for obtaining the value of a pseudopotential. It checks the type of the pseudopotential (`psp%type`).
    - If `psp%type == 1` (reciprocal space): It calls the internal function `PspRecipCubicSpline` to get the interpolated value.
    - If `psp%type == 2` (real space): This functionality is currently disabled and will result in an error.
  - **Arguments:**
    - `psp :: TYPE(pseudoPot), INTENT(IN)`: The pseudopotential object (defined in `CellInfo`) containing the tabulated data and spline coefficients.
    - `lookUp :: REAL(kind=DP), INTENT(IN)`: The magnitude of the q-vector (for reciprocal space) or r-vector (for real space, if enabled) at which to evaluate the pseudopotential.
  - **Return Value:**
    - `value :: REAL(kind=DP)`: The interpolated value of the pseudopotential at `lookUp`.

- **`FUNCTION PspRecipCubicSpline(psp, qNorm) RESULT(interpolated_value)`** (Internal to `PseudoPotLookup`):
  - **Description:** Performs cubic spline interpolation for a reciprocal space pseudopotential.
    1.  Handles edge cases: If `qNorm == 0.0`, it returns `psp%potValues(1)` (the first tabulated value, assumed to be the q=0 limit of the non-Coulombic part). If `qNorm >= psp%maxG` (beyond the tabulated range), it returns `0.0`.
    2.  For `0 < qNorm < psp%maxG`, it determines the position (`pos`) of `qNorm` within the tabulated `q` grid (derived from `psp%maxG` and `psp%numPoints`).
    3.  It then performs a cubic spline interpolation directly (without calling `spline_cubic_val` from `MathSplines` as an external function, but implementing its logic) using the pre-calculated second derivatives stored in `psp%potDD` and the tabulated values `psp%vqS` (which represent `V(q) + 4*pi*Z/q^2`).
    4.  Finally, it subtracts the Coulomb tail `4*pi*Z/qNorm^2` from the interpolated result to get the actual `V(q)`.
  - **Arguments:**
    - `psp :: TYPE(pseudoPot), INTENT(IN)`: The pseudopotential object.
    - `qNorm :: REAL(KIND=DP), INTENT(IN)`: The magnitude of the q-vector.
  - **Return Value:**
    - `interpolated_value :: REAL(KIND=DP)`: The value of `V(q)` at `qNorm`.

- **`FUNCTION PseudoPotDiffLookup(psp, qNorm) RESULT(derivative_value)`**:
  - **Description:** Calculates the first derivative `dV(q)/dq` of a reciprocal space pseudopotential at a given `qNorm`.
    1.  Checks if `qNorm > psp%maxG`, returning `0.0` if true.
    2.  Otherwise, it uses `spline_cubic_val` (from `MathSplines`) to find the derivative of the modified pseudopotential `psp%vqS` (i.e., `d/dq [V(q) + 4*pi*Z/q^2]`). The `spline_cubic_val` function returns this derivative as `ypval`.
    3.  It then adds the analytical derivative of the Coulomb tail `d/dq [-4*pi*Z/q^2] = +8*pi*Z/q^3` to `ypval / pspSpacing` (where `pspSpacing` is `dq`) to get the final `dV(q)/dq`.
  - **Arguments:**
    - `psp :: TYPE(pseudoPot), INTENT(IN)`: The pseudopotential object.
    - `qNorm :: REAL(KIND=DP), INTENT(IN)`: The magnitude of the q-vector.
  - **Return Value:**
    - `derivative_value :: REAL(KIND=DP)`: The value of `dV(q)/dq` at `qNorm`.

# Important Variables/Constants

The functionality of this module relies heavily on the structure and content of the `TYPE(pseudoPot)` object, which is defined in the `CellInfo` module. Key fields from `psp` used here include:
- `psp%type :: INTEGER`: Determines if the pseudopotential is real or reciprocal space.
- `psp%maxG :: REAL(KIND=DP)`: The maximum `q` value in the tabulated data for reciprocal space PSPs.
- `psp%numPoints :: INTEGER`: The number of data points in the tabulation.
- `psp%charge :: INTEGER`: The ionic charge `Z`, used for the Coulomb tail `4*pi*Z/q^2`.
- `psp%t(:) :: REAL(KIND=DP), POINTER`: Pointer to the array of `q` (or `r`) values (knots for the spline).
- `psp%vqS(:) :: REAL(KIND=DP), POINTER`: Pointer to the tabulated values of the modified pseudopotential `V(q) + 4*pi*Z/q^2` for reciprocal space PSPs.
- `psp%potDD(:) :: REAL(KIND=DP), POINTER`: Pointer to the pre-calculated second derivatives of `psp%vqS` at the knots, used for cubic spline interpolation.
- `psp%potValues(:) :: REAL(KIND=DP), POINTER`: Pointer to the original tabulated pseudopotential values. `psp%potValues(1)` is used as the `q=0` limit in `PspRecipCubicSpline`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter) and `PI`.
- **`MathSplines`**: For the `spline_cubic_val` function, which is used in `PseudoPotDiffLookup` to get the value and derivative from a cubic spline. While `PspRecipCubicSpline` implements its own spline evaluation, it relies on `psp%potDD` which would have been set up using `spline_cubic_set` from `MathSplines` when the pseudopotential was initially processed (likely in `ReadIonFile` or `ReadPseudo`).
- **`CellInfo`**: For the definition of `TYPE(pseudoPot)`.

This module is called by other parts of the code (e.g., `IonElectron` module when calculating `V_ei(G)`) whenever the value of a local pseudopotential or its derivative is needed at a specific `q` point that may not exactly coincide with the tabulated grid points.
