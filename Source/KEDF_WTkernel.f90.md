# Overview

The `KEDF_WTkernel` module is responsible for setting up and initializing the non-local kinetic energy kernel used primarily by the Wang-Teter (WT) KEDF, but the `keKernel` array it manages is also utilized by other non-local KEDFs like Wang-Govind-Carter (WGC) and Huang-Carter 2010 (HC10).

For the WT KEDF, the kernel is based on the Lindhard response function, `LindG(eta, lambda_param, mu_param)`, where `eta = |G| / (2 * k_F(rho0))`, and `k_F(rho0)` is the Fermi wavevector of a reference uniform electron gas density `rho0`. The parameters `lambda_param` and `mu_param` are typically the coefficients of the Thomas-Fermi and von Weizsacker functionals, respectively.

This module provides routines to:
1.  Compute this kernel in reciprocal space (`FillWT_ReciprocalSpace`).
2.  Read a pre-computed 1D kernel from a file and interpolate it onto the current reciprocal space grid (`FillWT_RestartFrom_GSpace`).

The computed kernel is stored in the first component (`keKernel(:,:,:,1)`) of the module-level `keKernel` array.

# Key Components

- **`MODULE KEDF_WTkernel`**: The main container for WT kernel generation.

- **`SUBROUTINE FillWT()`**:
  - **Description:** This is the primary public interface for generating the WT kernel. It first calculates two important parameters based on `rho0` (a reference density from `SYS` module) and the density exponents `alpha` and `beta`:
    - `coef = 5.0 / (9.0 * alpha * beta * rho0**(alpha + beta - 5.0/3.0))`
    - `tkF = 2.0 * (3.0 * rho0 * PI**2)**(1.0/3.0)` (this is `2 * k_F(rho0)`)
    It then calls `FillWT_ReciprocalSpace` to populate the kernel array.
  - **Key actions:** Calculates `coef` and `tkF`, then calls `FillWT_ReciprocalSpace`.

- **`SUBROUTINE FillWT_ReciprocalSpace()`**:
  - **Description:** This routine fills the `keKernel(:,:,:,1)` array with the Wang-Teter kernel values. For each reciprocal space vector `G` (represented by `qTable(ix,i2,i3)` which stores `|G|`), it calculates `eta = |G| / tkF` and then computes `keKernel(ix,i2,i3,1) = LindG(eta, lambda_tf, mu_vw) * coef`. (`lambda_tf` and `mu_vw` are imported from `KEDF_TF` and `KEDF_VW` respectively). It also includes logic to output the 1D kernel `LindG(eta, lambda_tf, mu_vw)` to a file (`KEDF-kernel-G.dat`) if the `outputKernel` flag (from `OUTPUT` module) is set.
  - **Output:** Populates `keKernel(:,:,:,1)`.

- **`SUBROUTINE FillWT_RestartFrom_GSpace(kernelFile, kernel, qNorm)`**:
  - **Description:** Reads a pre-computed 1D Wang-Teter kernel from a specified `kernelFile`. The file is expected to contain pairs of `eta_in_file` and `kernel_value_in_file`. This routine then interpolates these tabulated values (after scaling the kernel values by `coef`) onto the current simulation's reciprocal space grid. For each `eta = qNorm(ix,i2,i3) / tkF` of the current grid, it finds the corresponding kernel value by linear interpolation from the `eta_in_file`, `w0_in` data read from the file. The result is stored in `kernel(:,:,:,1)` (which is `keKernel(:,:,:,1)`).
  - **Arguments:**
    - `kernelFile :: CHARACTER(LEN=*), INTENT(IN)`: Name of the file containing the 1D kernel data.
    - `kernel :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`: The `keKernel` array to be filled.
    - `qNorm :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Magnitudes of G-vectors for the current grid (`qTable`).

# Important Variables/Constants

- **Module-Level Parameters (Must be set before calling `FillWT`):**
    - `alpha :: REAL(KIND=DP)`: Exponent for the `rho^alpha` term in WT/WGC functionals (default: -100.0, indicating it needs to be set, e.g., by `SetupKEDF`).
    - `beta :: REAL(KIND=DP)`: Exponent for the `rho^beta` term (default: -100.0, needs to be set).
- **Module-Level Kernel Storage:**
    - `keKernel :: REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE`: The array where the computed kernel is stored. `keKernel(:,:,:,1)` is used for the WT kernel. Other KEDF modules (like `KEDF_WGCkernel`, `KEDF_MGPkernel`) also populate other parts or re-purpose this array. This array must be allocated before `FillWT` is called (typically in `SetupKEDF`).
- **Calculated Parameters (in `FillWT`):**
    - `coef :: REAL(KIND=DP)`: A prefactor for the WT kernel, dependent on `alpha`, `beta`, `rho0`.
    - `tkF :: REAL(KIND=DP)`: `2 * k_F(rho0)`, where `k_F` is the Fermi wavevector of the reference density `rho0`.
- **Unused/Legacy Real-Space Kernel Variables:**
    - `rinc, maxr4K2, wtKernelr, wtKernelFr`: These appear to be related to an older, real-space method for kernel generation and are not used by the reciprocal space routines.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision) and `PI`.
- **`FOURIER`**: The `FFT` interface is imported but not directly used in the kernel generation routines themselves. It would be used by the KEDF modules (like `KEDF_WT`) that consume the generated kernel.
- **`SYS`**: For `rho0` (reference density for `tkF` and `coef`) and `rhoS` (reference density, though `rho0` is primarily used for WT kernel parameters here; `rhoS` is more relevant for WGC kernel generation which might also use/share `keKernel`).
- **`MathFunctions`**: For the `LindG` function (Lindhard response function), which forms the core of the WT kernel.
- **`PlaneWave`**: For `qTable` (magnitudes of G-vectors).
- **`MPI_Functions`**: For `rankGlobal` (to control file output by only one rank) and potentially utilities used by `WrtOut`.
- **`OutputFiles`, `OUTPUT`**: For `outputUnit`, `errorUnit`, `WrtOut`, and the `outputKernel` flag (controls diagnostic output of the 1D kernel).
- **`KEDF_TF`**: For `lambda` (TF coefficient, passed as a parameter to `LindG`).
- **`KEDF_VW`**: For `mu` (vW coefficient, passed as a parameter to `LindG`).

**Workflow:**
1.  The parameters `alpha`, `beta` (module level) and `rho0` (from `SYS`) must be set. `keKernel` must be allocated.
2.  `FillWT` is called. It calculates `coef` and `tkF`.
3.  `FillWT_ReciprocalSpace` is then called, which populates `keKernel(:,:,:,1)` by evaluating `LindG(qTable/tkF, lambda, mu) * coef` for each G-vector.
4.  Alternatively, `FillWT_RestartFrom_GSpace` can be used to load the kernel from a file.
5.  The `KEDF_WT` module (and potentially other KEDFs if they use the `keKernel(:,:,:,1)` component) then uses this populated `keKernel` for calculating energy and potential.

This module is central to defining the non-local behavior of the Wang-Teter KEDF and provides the necessary kernel data for it. The sharing of `keKernel` with other modules like `KEDF_WGCkernel` suggests a design where `keKernel` is a common storage for different types of non-local kernels.
