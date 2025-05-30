# Overview

The `KEDF_WGCkernel` module is responsible for generating and initializing the non-local kinetic energy kernel specifically for the Wang-Govind-Carter (WGC) KEDF. The WGC KEDF, detailed in Phys. Rev. B. 60, 16350 (1999), utilizes a density-dependent kernel to capture non-local effects in the kinetic energy of an electron system.

This module provides routines to:
1.  Compute the WGC kernel components by solving its defining Ordinary Differential Equation (ODE) using an external solver (via `IntKernelODE` module). The solution `w(eta)` and its relevant derivatives are then interpolated onto the reciprocal space grid.
2.  Alternatively, read a pre-computed WGC kernel from a file and interpolate it onto the current grid.

The computed kernel components are stored in the `keKernel` array, which is imported from the `KEDF_WTkernel` module, suggesting a shared storage or common interface for different non-local kernels.

# Key Components

- **`MODULE KEDF_WGCkernel`**: The main container for WGC kernel generation.

- **`SUBROUTINE FillWGC()`**:
  - **Description:** This is the primary public interface for generating the WGC kernel. It calculates `tkFstar = 2 * k_F(rhoS)` (where `rhoS` is a reference density) and then calls `FillWGC_ReciprocalSpace` to perform the actual computation and storage of the kernel.
  - **Key actions:** Sets `tkFstar`, calls `FillWGC_ReciprocalSpace`.

- **`SUBROUTINE FillWGC_ReciprocalSpace()`**:
  - **Description:** This core routine computes the WGC kernel components for each reciprocal space vector `G`.
    1.  If the kernel has not been computed before (`have_RK == .FALSE.`):
        a.  Sets parameters for the WGC ODE (`ode_alpha = alpha`, `ode_beta = beta`, `ode_gamma = gamma`, `wInf`, `KEDFtype = 2`).
        b.  Calls `makeKernel2` (from `IntKernelODE`) to solve the 2nd order ODE for `w(eta)`, `w'(eta)`, and `w''(eta)`. The results from `IntKernelODE` are `int_eta`, `int_w`, `int_w1` (for `eta*w'`), `int_w2` (for `eta^2*w''`).
        c.  Sets up spline coefficients (`nls_wpp`, `nls_w1pp`, `nls_w2pp`) for these solutions using `spline_cubic_set`.
        d.  Sets `have_RK = .TRUE.`.
    2.  Iterates over all G-vectors (`qTable(ix,i2,i3)`):
        a.  Calculates `eta = |G| / tkFstar`.
        b.  Uses `spline_cubic_val` to interpolate `int_w`, `int_w1`, `int_w2` at the current `eta` to get `keKernel(ix,i2,i3,1)` (for `w(eta)`), `keKernel(ix,i2,i3,2)` (for `eta*w'(eta)`), and `keKernel(ix,i2,i3,3)` (for `eta^2*w''(eta)`).
        c.  Derives `keKernel(ix,i2,i3,4)` from these components and `gamma`, `rhoS`.
        d.  Transforms `keKernel(ix,i2,i3,2)` and `keKernel(ix,i2,i3,3)` to store different forms of derivatives.
    3.  Applies post-processing to `keKernel` components based on the `WGCT` flag (e.g., zeroing out higher-order terms if `WGCT` indicates a simpler version of the WGC functional).
  - **Output:** Populates the `keKernel` array with WGC-specific values.

- **`SUBROUTINE FillWGC_RestartFrom_ReciprocalSpace(kernelFile, kernel, qNorm)`**:
  - **Description:** Reads a pre-computed WGC kernel from a specified `kernelFile`. The file is expected to contain `eta`, `w(eta)`, `w'(eta)`, and `w''(eta)`. This routine then interpolates these tabulated values onto the current simulation's reciprocal space grid (`qNorm`, which is `qTable`) using splines (implicitly, as `nls_wpp` etc. must be set up before calling `spline_cubic_val` if this routine were to use them; however, it seems to directly use the read-in arrays `eta_in`, `w0_in`, etc. for interpolation with `spline_cubic_val` assuming `nls_wpp` etc. are already available from a previous `FillWGC_ReciprocalSpace` call if `have_RK` is true). It then processes these interpolated values similar to `FillWGC_ReciprocalSpace` to populate the `kernel` array (which is `keKernel`).
  - **Arguments:**
    - `kernelFile :: CHARACTER(LEN=*), INTENT(IN)`: Name of the file containing the kernel data.
    - `kernel :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`: The `keKernel` array to be filled.
    - `qNorm :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Magnitudes of G-vectors for the current grid (`qTable`).

# Important Variables/Constants

- **Module-Level Parameters for WGC Kernel:**
    - `gamma :: REAL(KIND=DP)`: An averaging parameter for the WGC functional (default: -1.0). Used in deriving `keKernel(:,:,:,4)`.
    - `firstOrderWGC :: INTEGER`: Flag indicating the order of Taylor expansion for WGC (default: -1, which usually implies a second-order or more complete form via `WGCT=2`). Not directly used in kernel generation but by `KEDF_WGC`.
    - `WGCT :: INTEGER`: A flag that controls which terms of the WGC kernel are used or how they are processed (default: 2). Affects post-processing of `keKernel`.
    - `alpha5, beta5 :: REAL(KIND=DP)`: Coefficients used if `WGCT = -5` (default: -100.0, effectively unset).
    - `have_RK :: LOGICAL`: Flag, set to `.TRUE.` after the WGC kernel's ODE has been solved for the first time. Prevents re-solving the ODE unnecessarily.
- **Spline Coefficients (Module-Level, Allocatable):**
    - `nls_wpp(:), nls_w1pp(:), nls_w2pp(:) :: REAL(KIND=DP), ALLOCATABLE`: Store second derivatives for spline interpolation of `int_w`, `int_w1`, and `int_w2` respectively.
- **Imported Variables:**
    - `keKernel :: REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:,:)` (from `KEDF_WTkernel`): The array where this module stores the computed WGC kernel components.
    - `rhoS :: REAL(KIND=DP)` (from `SYS`): Reference density used for calculating `tkFstar` and in kernel component formulas.
    - `alpha, beta :: REAL(KIND=DP)` (from `KEDF_WTkernel`): Density exponents, passed to `IntKernelODE` as `ode_alpha`, `ode_beta`.
    - `qTable :: REAL(KIND=DP), DIMENSION(:,:,:)` (from `PlaneWave`): Magnitudes of G-vectors.
    - `outputKernel :: LOGICAL` (from `OUTPUT`): Flag to control printing of kernel data to a file.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`, `PI`.
- **`FOURIER`**: `FFT` is imported but not directly used in these routines.
- **`PlaneWave`**: For `qTable`.
- **`KEDF_WTkernel`**: Crucially for the `keKernel` array (where the output is stored) and for `alpha`, `beta` parameters. This indicates a shared infrastructure for different non-local KEDFs.
- **`SYS`**: For `rho0` (though `rhoS` is directly used) and `rhoS`.
- **`MPI_Functions`**: For `rankGlobal` (for rank-specific output) and potentially utilities used by `WrtOut`.
- **`OutputFiles`, `OUTPUT`**: For `outputUnit`, `WrtOut`, `message`, and the `outputKernel` flag.
- **`IntKernelODE`**: This is a primary dependency for `FillWGC_ReciprocalSpace`. It uses:
    - Control parameters: `KEDFType`, `winf`, `ode_alpha`, `ode_beta`, `ode_gamma`.
    - ODE solution results: `int_eta` (eta grid), `int_w` (w(eta)), `int_w1` (eta*w'), `int_w2` (eta^2*w'').
    - The ODE solver driver: `makeKernel2`.
    - `IntKernelODEClean` (aliased as `Clean`, though not called directly in these routines, presumably called elsewhere).
- **`MathSplines`**: For `spline_cubic_set` (to prepare spline coefficients) and `spline_cubic_val` (to interpolate kernel values).

**Workflow:**
1.  `FillWGC` is called by a higher-level KEDF setup routine if WGC is selected.
2.  `FillWGC_ReciprocalSpace` is invoked.
3.  If `have_RK` is false (first time or after a system change requiring re-computation):
    a.  It sets parameters for the WGC ODE in `IntKernelODE`.
    b.  Calls `makeKernel2` from `IntKernelODE`, which solves the ODE and stores the solution (`int_w`, `int_w1`, `int_w2` on an `int_eta` grid).
    c.  Sets up spline coefficients (`nls_wpp`, etc.) for these solutions.
    d.  Sets `have_RK = .TRUE.`.
4.  It then iterates through all G-points, calculates `eta = |G| / tkFstar`, and uses `spline_cubic_val` to interpolate the ODE solutions to get kernel values at this `eta`.
5.  These interpolated values are processed and stored in the `keKernel` array.
6.  `keKernel` is post-processed based on the `WGCT` flag.
Alternatively, if restarting, `FillWGC_RestartFrom_ReciprocalSpace` can be called to load and interpolate kernel data from a file. The `KEDF_WGC` module then uses this `keKernel` for its calculations.
