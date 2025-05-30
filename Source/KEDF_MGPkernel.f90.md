# Overview

The `KEDF_MGPkernel` module is dedicated to the construction of the nonlocal kinetic energy kernel for the Mi-Genova-Pavanello (MGP) Kinetic Energy Density Functional (KEDF). The MGP KEDF is detailed in J. Chem. Phys. 148, 184107 (2018).

This module calculates the kernel components in reciprocal space (`G-space`) and stores them in the `keKernel` array, which is then used by the `KEDF_MGP` module to compute the actual MGP kinetic energy and potential.

# Key Components

- **`MODULE KEDF_MGPkernel`**: The main container for the MGP kernel generation routines.

- **`SUBROUTINE FillMGP()`**:
  - **Description:** This is the primary public interface for generating the MGP kernel. It initializes a couple of constants (`cMGP`, `tkF` based on `rho0`) and then calls `FillMGP_ReciprocalSpace` to perform the actual kernel calculation.
  - **Key actions:** Sets `cMGP` and `tkF`, then calls `FillMGP_ReciprocalSpace`.

- **`SUBROUTINE FillMGP_ReciprocalSpace()`**:
  - **Description:** This core subroutine computes the different components of the MGP kernel for each reciprocal space vector `G` (represented by `qTable(ix,i2,i3)`). The MGP kernel in this implementation appears to have three main parts stored in different indices of the `keKernel` array:
    1.  **`keKernel(ix,i2,i3,1)` (Total Kernel):** This is the sum of a term proportional to the Lindhard function `LindG(q/tkF,1.0,1.0)` and the other two calculated components (`keKernel(:,:,:,2)` and `keKernel(:,:,:,3)`).
    2.  **`keKernel(ix,i2,i3,2)` (Hypercorrelation Term):** Calculated via a numerical integration involving the Lindhard function `LindG` with specific parameters (-0.6, 1.0). The integration is over a variable `t_var`.
    3.  **`keKernel(ix,i2,i3,3)` (Kinetic Electron Term):** This term has a functional form involving `erf(q)^2 * LumgpFactor * exp(-q^2*LumgpEXP)/q^2`.
    The final kernel components are scaled by `cMGP`. The `q=0` limit is handled by setting the kernel values to zero.
  - **Output:** Populates the `keKernel` array (imported from `KEDF_WTkernel`).

# Important Variables/Constants

- **Module-Level Variables (calculated in `FillMGP`):**
    - `tkF :: REAL(KIND=DP)`: A wavevector scale, calculated as `2.0 * (3.0 * rho0 * PI**2)**(1.0/3.0)`.
    - `cMGP :: REAL(KIND=DP)`: A constant, `PI**2 / (3.0 * PI**2)**(1.0/3.0)`.

- **Imported Variables Used in Kernel Calculation:**
    - `rho0 :: REAL(KIND=DP)` (from `SYS`): A reference electron density used to calculate `tkF`.
    - `LumgpExp :: REAL(KIND=DP)` (from `SYS`): An exponent factor used in the "Kinetic Electron Term."
    - `LumgpFactor :: REAL(KIND=DP)` (from `SYS`): A multiplicative factor used in the "Kinetic Electron Term."
    - `qTable :: REAL(KIND=DP), DIMENSION(:,:,:)` (from `PlaneWave`): Array containing the magnitudes of the reciprocal space vectors `q = |G|`.
    - `keKernel :: REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:,:)` (from `KEDF_WTkernel`): The module-level array where the calculated MGP kernel components are stored. This array is expected to be allocated before `FillMGP` is called.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision) and `PI`.
- **`FOURIER`**: The `FFT` interface is imported but not explicitly used within the `FillMGP_ReciprocalSpace` routine.
- **`SYS`**: For system parameters `rho0`, `LumgpExp`, and `LumgpFactor`.
- **`MathFunctions`**: For the `LindG` function (Lindhard response function).
- **`PlaneWave`**: For `qTable` (magnitudes of G-vectors).
- **`MPI_Functions`**: Module is used, but no direct MPI calls are visible in these routines. Utilities from it might be used by `WrtOut`.
- **`OutputFiles`**: Module is used (likely for `outputUnit` used by `WrtOut`).
- **`OUTPUT`**: For `WrtOut` (writing messages) and `outputUnit`.
- **`KEDF_WTkernel`**: This is a crucial interaction, as `KEDF_MGPkernel` fills the `keKernel` array that is declared and likely used by `KEDF_WTkernel` or other KEDF modules that might share/store kernels there. This suggests a design where `KEDF_WTkernel` might be a more general module for handling various types of kernels.
- **`KEDF_TF`**: Imports `lambda`, but it's not used.
- **`KEDF_VW`**: Imports `mu`, but it's not used.
- **`CellInfo`**: Imports `cell`, but it's not used.

**Workflow:**
1.  The `FillMGP` subroutine is called, typically during the setup phase of a KEDF calculation if the MGP functional is selected.
2.  `FillMGP` initializes `tkF` and `cMGP`.
3.  `FillMGP_ReciprocalSpace` is then called.
4.  This routine iterates over all G-vectors (grid points in reciprocal space).
5.  For each G-vector, it computes the three components of the MGP kernel using the formulas involving `qTable`, `tkF`, `LindG`, `LumgpExp`, and `LumgpFactor`.
6.  The results are stored in the slices `keKernel(:,:,:,1)`, `keKernel(:,:,:,2)`, and `keKernel(:,:,:,3)`.
7.  The `KEDF_MGP` module then uses `keKernel(:,:,:,1)` (the total kernel) for its energy and potential calculations.

The module effectively populates a pre-defined array (`keKernel`) with the MGP-specific kernel values, which are then consumed by the `KEDF_MGP` module.
