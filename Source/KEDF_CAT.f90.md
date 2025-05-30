# Overview

The `KEDF_CAT` module implements the Corrected Average Thomas-Fermi (CAT) Kinetic Energy Density Functional (KEDF). This functional is based on the work of David Garcia-Aldea and J.E. Alvarellos, as detailed in Phys. Rev. A 76, 052504 (2007). The CAT KEDF aims to provide a more accurate description of the kinetic energy of an electron system compared to simpler local or semi-local approximations, by incorporating non-local information through a kernel.

The module provides routines to:
- Calculate the CAT kinetic energy and its corresponding potential.
- Include vacuum corrections for systems with regions of very low electron density.
- Prepare the necessary non-local kernel by solving a differential equation and splining the solution.

# Key Components

- **`MODULE KEDF_CAT`**: The main container for CAT KEDF functionalities.

- **`SUBROUTINE FillCAT`**:
  - **Description:** This crucial subroutine prepares the non-local kinetic energy kernel (`keKernel`) required by the CAT functional. It sets parameters (`ode_alpha`, `ode_beta`, `ode_gamma` which are linked to `cat_alpha`, `cat_beta`, `cat_gamma`) for an ordinary differential equation (ODE), solves this ODE using routines from the `IntKernelODE` module (specifically `makeKernel2`), and then uses splines (`MathSplines` module) to interpolate the solution (`int_w`, `int_w1`, `int_w2`) onto the reciprocal space grid points, storing different components or derivatives in `keKernel`.
  - **Key actions:** Sets `KEDFType=1`, defines CAT parameters, calls `makeKernel2`, and then populates `keKernel` using spline interpolation of the ODE solution.

- **`SUBROUTINE CalCAT(potential, rho, calcEnergy, optSqrt, bvac, eTF, eVW, eNL)`**:
  - **Description:** This is the main driver routine for calculating the CAT KEDF potential and its energy components. It combines a scaled Thomas-Fermi (TF) term, the von Weizsaecker (vW) term, and the non-local CAT term (which itself might be a correction or addition to a base TF-like term). It can optionally apply vacuum corrections.
  - **Arguments:**
    - `potential :: REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(OUT)`: Output KEDF potential.
    - `rho :: REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(IN)`: Input electron density.
    - `calcEnergy :: LOGICAL, INTENT(IN)`: Flag to calculate energy components.
    - `optSqrt :: LOGICAL, INTENT(IN)`: Flag indicating if optimization is w.r.t. `sqrt(rho)`.
    - `bvac :: LOGICAL, INTENT(IN)`: Flag to apply vacuum corrections.
    - `eTF, eVW, eNL :: REAL(KIND=DP), INTENT(OUT)`: Output energy components (TF, vW, Non-Local/CAT).

- **`FUNCTION CATEnergy(rho) RESULT(energy_val)`**:
  - **Description:** Calculates the non-local part of the CAT kinetic energy, using `CATrhot` to get an intermediate density-dependent quantity.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION CATPot(rho) RESULT(potential_val)`**:
  - **Description:** Calculates the potential corresponding to `CATEnergy`. Involves FFTs of various density-derived terms convolved with components of `keKernel`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION CATrhot(rho) RESULT(rho_tilde)`**:
  - **Description:** Calculates an intermediate quantity `rho_tilde` based on `rho` and convolutions with `keKernel`. This `rho_tilde` is then used in `CATEnergy` and `CATPot`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION CATvacEnergy(rho) RESULT(energy_val)`**:
  - **Description:** Similar to `CATEnergy` but incorporates vacuum corrections using `CutoffFunc`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION CATvacPot(rho) RESULT(potential_val)`**:
  - **Description:** Similar to `CATPot` but incorporates vacuum corrections using `CutoffFunc` and `CutoffFuncD`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION CATvacRhot(rho) RESULT(rho_tilde_vac)`**:
  - **Description:** Similar to `CATrhot` but incorporates vacuum corrections.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`SUBROUTINE CATPotentialPlus(rho, potential, calcEnergy, energy, vacCutoff)`**:
  - **Description:** An alternative routine for calculating the CAT potential and energy, mentioned in comments as "still being debugged."
  - **Arguments:** Similar to `CalCAT` but with slightly different output structure.

# Important Variables/Constants

- **Module-Level Parameters for CAT KEDF:**
    - `cat_alpha :: REAL(KIND=DP)`: Default 0.0.
    - `cat_beta :: REAL(KIND=DP)`: Default 2.0/3.0.
    - `cat_gamma :: REAL(KIND=DP)`: Default 1.4. These parameters define the specific form of the CAT functional and are used to set `ode_alpha`, `ode_beta`, `ode_gamma` for the kernel ODE solver.
- **Shared Variables from `KEDF_WGC` (via `USE KEDF_WGC`):**
    - `rhoS :: REAL(KIND=DP)`: A reference density used in kernel generation.
    - `cTF :: REAL(KIND=DP)`: The Thomas-Fermi constant `(3/10)*(3*pi^2)^(2/3)`.
    - `keKernel :: REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:,:)`: The non-local kinetic energy kernel in reciprocal space, populated by `FillCAT`. Different indices of the 4th dimension store `w(eta)` and its scaled derivatives.
    - `CutoffFunc, CutoffFuncD`: Functions used for applying vacuum corrections (likely defined in `KEDF_WGC`).
    - `WGCT :: INTEGER`: An integer flag, likely from `KEDF_WGC`, that might control aspects of the kernel generation or its off-diagonal terms.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP`, `PI`.
- **`OutputFiles`**: For `outputUnit`.
- **`KEDF_WGC`**: For shared variables like `rhoS`, `cTF`, `keKernel`, `CutoffFunc`, `CutoffFuncD`, and potentially control flags like `WGCT`. This indicates a close relationship or a hierarchical structure where CAT might be built upon or share components with WGC.
- **`Output`**: For `WrtOut`.
- **`KEDF_TF`**: For `TFPotential`, `TFEnergy` (to include the Thomas-Fermi part), and `lambda` (TF coefficient).
- **`KEDF_VW`**: For `VWPotentialSqrtPlus` (for the von Weizsaecker part), and `mu` (vW coefficient).
- **`CellInfo`**: For FFT grid dimensions (`n1G`, `n2G`, `n3G`, `k1G`, `k2G`, `k3G`) and the `cell` variable.
- **`Fourier`**: For the `FFT` interface, used extensively in `CATPot`, `CATrhot`, and their vacuum-corrected counterparts to perform convolutions via reciprocal space multiplication.
- **`PlaneWave`**: For `qTable` (G-vector magnitudes), used in `FillCAT` for kernel generation.
- **`IntKernelODE`**: This is a critical dependency for `FillCAT`. It provides:
    - Control parameters: `KEDFType`, `ode_alpha`, `ode_beta`, `ode_gamma`.
    - ODE solution results: `winf`, `int_eta` (eta grid), `int_w` (kernel function w(eta)), `int_w1` (w'), `int_w2` (w'').
    - The ODE solver driver: `makeKernel2`.
    - Cleanup routine: `IntKernelODEClean` (aliased as `Clean`).
- **`MathSplines`**: For `spline_cubic_set` and `spline_cubic_val`, used in `FillCAT` to interpolate the numerically solved kernel onto the simulation's G-grid.

The `KEDF_CAT` module relies on `FillCAT` to first generate its specific non-local kernel. Then, routines like `CalCAT` use this kernel along with TF and vW components to compute the total KEDF energy and potential. The vacuum-corrected functions provide modified behavior for systems with significant low-density regions.
