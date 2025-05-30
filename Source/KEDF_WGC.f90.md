# Overview

The `KEDF_WGC` module implements the Wang-Govind-Carter (WGC) nonlocal Kinetic Energy Density Functional (KEDF). This functional, detailed in Phys. Rev. B. 60, 16350 (1999), is a significant KEDF that uses a density-dependent kernel to improve upon local and semi-local approximations. It involves convolutions of powers of the electron density (`rho^alpha`, `rho^beta`) with components of this kernel.

The module provides routines to calculate the WGC kinetic energy potential (`WGCPotentialPlus`), the WGC contribution to the stress tensor (`WGCStress`), and utility functions for applying a cutoff in vacuum regions (`CutoffFunc`, `CutoffFuncD`). The specific form of the WGC functional (e.g., first or second order Taylor expansion) is controlled by flags imported from `KEDF_WGCkernel`.

# Key Components

- **`MODULE KEDF_WGC`**: The main container for WGC KEDF functionalities.

- **`SUBROUTINE WGCPotentialPlus(rho, potential, calcEnergy, energy, vacCutoff)`**:
  - **Description:** Calculates the WGC potential and, optionally, the energy. The calculation depends on the `firstOrderWGC` flag (for 1st or 2nd order Taylor expansion of the WGC functional) and potentially the `WGCT` flag (for different kernel term combinations). It involves FFTs of `rho^alpha` and `rho^beta` (or `rho^alpha*(rho-rhoS)`, etc.) convolved with different components of the `keKernel`. If `vacCutoff` is true, `CutoffFunc` and its derivative `CutoffFuncD` are applied to density terms, particularly `rho^beta`.
  - **Arguments:**
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density.
    - `potential :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output WGC potential contribution.
    - `calcEnergy :: LOGICAL, INTENT(IN)`: If `.TRUE.`, calculates the WGC energy.
    - `energy :: REAL(KIND=DP), INTENT(OUT)`: Calculated WGC energy.
    - `vacCutoff :: LOGICAL, INTENT(IN)`: If `.TRUE.`, applies vacuum corrections.

- **`FUNCTION CutoffFunc(rho) RESULT(cutoff_val)`**:
  - **Description:** A cutoff function applied in vacuum regions to ensure numerically stable behavior of the functional where density is very low. The form is `(exp(rho/rhoStep) - 1) / (exp(rho/rhoStep) + exp(rhoV/rhoStep))`, with Taylor expansions for small `rho/rhoStep`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.
  - **Return Value:** `cutoff_val :: REAL(KIND=DP), DIMENSION(n1G, n2G, n3G)`.

- **`FUNCTION CutoffFuncD(rho) RESULT(deriv_cutoff_val)`**:
  - **Description:** Calculates the derivative of `CutoffFunc` with respect to `rho`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION CutoffFuncOld(rho) RESULT(cutoff_val_old)` / `FUNCTION CutoffFuncDold(rho) RESULT(deriv_cutoff_val_old)`**:
  - **Description:** Older versions of the cutoff function and its derivative. `CutoffFunc` includes a bug fix for the Taylor expansion compared to `CutoffFuncOld`.

- **`FUNCTION WGCStress(rhoReal) RESULT(stress_tensor)`**:
  - **Description:** Calculates the WGC KEDF contribution to the stress tensor. This is a complex calculation involving FFTs of various powers of the density (`rho^alpha`, `rho^beta`, `rho^(alpha+1)`, etc.) convolved with different kernel components and their derivatives (implicitly, via terms like `c00`, `c01`, etc., which depend on kernel components and other parameters like `alpha`, `beta`, `gamma`, `rhoS`, and flags `hold0`, `holdS`).
  - **Arguments:** `rhoReal :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

# Important Variables/Constants

- **Module-Level Parameters for Vacuum Cutoff:**
    - `rhoV :: REAL(KIND=DP)`: A small density value defining the vacuum region threshold (default: 1.0E-6).
    - `rhoStep :: REAL(KIND=DP)`: A density step parameter used in the cutoff function's exponential terms (default: 1.0E-6).
- **Imported Parameters and Kernel:**
    - `cTF :: REAL(KIND=DP)` (from `KEDF_TF`): The Thomas-Fermi constant.
    - `alpha, beta :: REAL(KIND=DP)` (from `KEDF_WTkernel`): Exponents for density powers in the WGC functional.
    - `keKernel :: REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:,:)` (from `KEDF_WTkernel`): The nonlocal kinetic energy kernel in reciprocal space. Different components (e.g., `keKernel(:,:,:,1)`, `keKernel(:,:,:,2)`, etc.) are used depending on the WGC formulation. This kernel is filled by routines in `KEDF_WGCkernel`.
    - `WGCT :: INTEGER` (from `KEDF_WGCkernel`): Flag that likely selects different terms or forms of the WGC kernel/functional.
    - `firstOrderWGC :: INTEGER` (from `KEDF_WGCkernel`): Flag to switch between first-order (`1`) and second-order (`-1`) Taylor-expanded forms of the WGC functional.
    - `gamma :: REAL(KIND=DP)` (from `KEDF_WGCkernel`): A parameter for the WGC kernel.
    - `rho0, rhoS :: REAL(KIND=DP)` (from `SYS`): Reference densities used in kernel generation and stress calculations.
    - `hold0, holdS :: LOGICAL` (from `SYS`): Flags indicating if `rho0` and `rhoS` are treated as constants with respect to cell volume changes (affects stress calculation).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP`.
- **`OutputFiles`**: For `outputUnit`.
- **`KEDF_TF`**: For `cTF`.
- **`KEDF_WTkernel`**: For `alpha`, `beta`, and the crucial `keKernel` array. This highlights a modular design where `KEDF_WGCkernel` prepares the kernel, which is then stored in a common place (accessible via `KEDF_WTkernel`) for `KEDF_WGC` to use.
- **`SetupFFT`**: Module is used, but no specific entities are called directly in the provided code; likely sets up global FFT parameters used by `Fourier_NEW`.
- **`Fourier_NEW`**: For `FFT_NEW` and `FFT_STD_STATE`, used for performing Fast Fourier Transforms.
- **`Output`**: For `WrtOut` and `outputUnit` (though `outputUnit` is also directly from `OutputFiles`).
- **`Fourier`**: The older `FFT` interface is imported but seems to be superseded by `Fourier_NEW` calls.
- **`PlaneWave`**: For `qTable` (G^2 values), `qMask`, and `qVectors`.
- **`CellInfo`**: For grid dimensions (`n1G`, `n2G`, `n3G`, `k1G`, `k2G`, `k3G`).
- **`KEDF_WGCkernel`**: For parameters `WGCT`, `firstOrderWGC`, and `gamma` that dictate the specific WGC functional form and kernel to be used. The actual kernel data is expected to be filled by this module into `keKernel`.
- **`SYS`**: For reference densities `rho0`, `rhoS` and flags `hold0`, `holdS`.
- **`MathFunctions`**: For `GPrime` (derivative of Lindhard function, used in `WGCStress`).

The `KEDF_WGC` module calculates one of the more complex, nonlocal KEDFs. Its operation relies significantly on:
1.  A pre-computed `keKernel` provided via `KEDF_WGCkernel` (and stored in `KEDF_WTkernel`).
2.  Control flags (`firstOrderWGC`, `WGCT`) also from `KEDF_WGCkernel`.
3.  FFT routines for convolutions.
4.  Optional vacuum correction functions.
