# Overview

The `KEDF_VW` module is responsible for calculating the von Weizsacker (vW) component of the kinetic energy, its corresponding potential, and its contribution to the stress tensor. The von Weizsacker KEDF is a gradient-dependent functional, considered exact for a system with only one or two electrons, and provides a lower bound for the kinetic energy of any system.

The vW kinetic energy is commonly expressed in two forms:
1.  In terms of the density `rho(r)`: `E_vW = mu * (1/8) * integral( |nabla rho(r)|^2 / rho(r) dr )`
2.  In terms of `phi(r) = sqrt(rho(r))`: `E_vW = mu * (1/2) * integral( |nabla phi(r)|^2 dr )`
    In reciprocal space, the second form becomes `E_vW = mu * (1/2) * Omega * SUM_G { G^2 * |phi(G)|^2 }`, where `phi(G)` is the Fourier transform of `phi(r)`.

The module provides functions based on these formulations, primarily using the `sqrt(rho)` representation for potential and energy calculations.

# Key Components

- **`MODULE KEDF_VW`**: The main container for vW KEDF calculations.

- **`SUBROUTINE CalVW(potential, rho, calcEnergy, optSqrt, eVW)`**:
  - **Description:** Calculates the vW potential and, optionally, the energy. It handles both spin-unpolarized and spin-polarized cases. For spin-polarized systems, it typically applies the vW functional to `2*rho_spin` for each spin channel and sums the energies (multiplied by 0.5). It uses `VWPotentialSqrtPlus` internally.
  - **Arguments:**
    - `potential :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT)`: Output vW potential (added to the input `potential`).
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: Input electron density (4th dimension for spin).
    - `calcEnergy :: LOGICAL, INTENT(IN)`: If `.TRUE.`, calculates the vW energy.
    - `optSqrt :: LOGICAL, INTENT(IN)`: If `.TRUE.`, the output `potential` is `dE/d(sqrt(rho))`; otherwise, it's `dE/drho`.
    - `eVW :: REAL(KIND=DP), INTENT(OUT)`: Calculated vW energy.

- **`FUNCTION VWEnergy(sqrtRhoR_SI) RESULT(energy_val)`**:
  - **Description:** Calculates the vW kinetic energy using the `sqrt(rho)` formulation in reciprocal space: `E_vW = 0.5 * mu * SUM(sqrtRhoR_SI * FFT_backward[ FFT_forward[sqrtRhoR_SI] * qTable^2 ])`. `qTable` here stores `G^2`.
  - **Arguments:** `sqrtRhoR_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Square root of spin-independent electron density in real space.

- **`FUNCTION VWPotential(rhoR_SI) RESULT(potential_val)`**:
  - **Description:** Calculates the vW potential `dE_vW/drho`. The implementation is effectively `-mu * (1/2) * laplacian(sqrt(rho)) / sqrt(rho)`. It computes this as `mu * 0.5 * FFT_backward[ FFT_forward[SQRT(rhoR_SI)] * qTable^2 ] / SQRT(rhoR_SI)`.
  - **Arguments:** `rhoR_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density.

- **`FUNCTION VWPotentialSqrt(rhoR_SI) RESULT(potential_sqrt_val)`**:
  - **Description:** Calculates the vW potential derivative with respect to `sqrt(rho)`, which is `dE_vW / d(sqrt(rhoR_SI)) = -mu * laplacian(sqrt(rhoR_SI))`. This is computed as `mu * FFT_backward[ FFT_forward[SQRT(rhoR_SI)] * qTable^2 ]`.
  - **Arguments:** `rhoR_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density.

- **`SUBROUTINE VWPotentialSqrtPlus(sqrtRho_SI, potentialSqrt, calcEnergy, energy, weightingfunc, dV, IfCalForce)`**:
  - **Description:** This is the primary routine used by `CalVW` and other KEDF modules (like `KEDF_Q`, `KEDF_EvW`) to get the vW potential (w.r.t `sqrt(rho)`) and energy. It calculates `potentialSqrt = mu * FFT_backward[ FFT_forward[sqrtRho_SI] * qTable^2 ]` and `energy = 0.5 * SUM(sqrtRho_SI * potentialSqrt) * mu` (before `potentialSqrt` is scaled by `mu`). Optional arguments for weighting functions are present but not used in the main path shown.
  - **Arguments:**
    - `sqrtRho_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Square root of spin-independent density.
    - `potentialSqrt :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output vW potential w.r.t. `sqrt(rho)`.
    - `calcEnergy :: LOGICAL, INTENT(IN)`: Flag to calculate energy.
    - `energy :: REAL(KIND=DP), INTENT(OUT)`: Output vW energy.
    - `weightingfunc, dV, IfCalForce`: Optional arguments, not fully utilized in the main logic.

- **`FUNCTION VWStress(sqrtRhoQ, sizeRealRho) RESULT(stress_tensor)`**:
  - **Description:** Calculates the vW contribution to the stress tensor. The formula used is `Stress_ab = -mu * SUM_G [ (G_a * FFT[sqrt(rho)]^*) * (G_b * FFT[sqrt(rho)]) ] / N_grid_real`.
  - **Arguments:**
    - `sqrtRhoQ :: COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Fourier transform of `sqrt(rho)`.
    - `sizeRealRho :: INTEGER, INTENT(IN)`: Total number of real-space grid points.
  - **Return Value:** `stress_tensor :: REAL(KIND=DP), DIMENSION(3,3)`.

# Important Variables/Constants

- **Module-Level Parameter:**
    - **`mu :: REAL(KIND=DP)`**: A scaling factor for the vW contribution. Default value is -100.0, which is unphysical. This parameter **must be initialized** to a physically meaningful value (e.g., 1.0 for full vW, or other fractions like 1/8, 1/9, 1/5 when vW is part of a composite KEDF) elsewhere in the code, likely during KEDF setup.
- **Key Reciprocal Space Quantities (from `PlaneWave` module):**
    - `qTable :: REAL(KIND=DP), DIMENSION(k1G,k2G,k3G)`: Array containing `G^2` (squared magnitudes of reciprocal lattice vectors).
    - `qVectors :: REAL(KIND=DP), DIMENSION(k1G,k2G,k3G,3)`: Array containing Cartesian components `(Gx, Gy, Gz)` of G-vectors.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision) and `IMAG` (imaginary unit).
- **`Fourier`**: For the older `FFT` interface (used in `VWPotential`, `VWEnergy`, `VWPotentialSqrt`).
- **`Fourier_NEW`**: For the newer `FFT_NEW` interface and `FFT_STD_STATE` (used in the primary workhorse routine `VWPotentialSqrtPlus`).
- **`PlaneWave`**: For `qTable` (G^2 values) and `qVectors` (G-vector components).
- **`OutputFiles`**: For `outputUnit` and `errorUnit` (Fortran I/O units).
- **`CellInfo`**: For `cell` (cell information), `numSpin` (used in `CalVW`), and grid dimensions (`n1G`, `n2G`, `n3G`, `k1G`, `k2G`, `k3G`).
- **`MathFunctions`**: `GlobalVal3d` is imported in `VWPotentialSqrtPlus` but not visibly used in the provided snippet.

The vW KEDF is a fundamental gradient correction. Routines from this module, especially `VWPotentialSqrtPlus`, are frequently called by other KEDF modules (e.g., `KEDF_Q`, `KEDF_EvW`, `KEDF_GGA` when `CP` is non-zero) to incorporate the vW energy and potential, often scaled by the `mu` parameter. The choice of FFT interface (`Fourier` vs. `Fourier_NEW`) varies between older functions and the more recent `VWPotentialSqrtPlus`.
