# Overview

The `XC_PBE` module implements the Perdew-Burke-Ernzerhof (PBE) Generalized Gradient Approximation (GGA) for the exchange-correlation (XC) energy, potential, and stress tensor. PBE is a widely used functional in Density Functional Theory (DFT) calculations.

This module provides two main ways to calculate PBE quantities:
1.  Direct implementation of the PBE formulas (e.g., in `PBEPot` for spin-unpolarized cases and `PBEStress`).
2.  An interface to the LibXC library (`PBE_LibXC`) for a more general handling of PBE, including spin-polarized cases. LibXC is an external library of exchange-correlation functionals.

**References:**
- J.P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996) (The original PBE paper).
- LibXC library documentation (for `PBE_LibXC` functionality).

# Key Components

- **`MODULE XC_PBE`**: The main container for PBE XC calculations.

- **Module-Level Parameter:**
    - `pbeCutoff :: REAL(KIND=DP)`: A cutoff value for the squared gradient of density (default: 0.0). If `|nabla rho|^2` is less than this cutoff, the gradient is treated as zero in the `PBEPot` subroutine. This is intended to improve numerical stability in very low-density or slowly varying regions, though the comment notes it "should not use it in principle."

- **PBE Calculation Routines:**
    - **`SUBROUTINE PBEPot(rhoReal_SI, potential, calcEnergy, energy)`**:
      - **Description:** Calculates the PBE XC potential and optionally the energy for **spin-unpolarized** systems using a direct implementation of the PBE formulas. It involves calculating the reduced density gradient `s`, the PBE exchange enhancement factor `Fx(s)`, and the PBE correlation energy density `ec(rs, zeta=0, t)` and their derivatives.
      - **Arguments:**
        - `rhoReal_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent real-space electron density.
        - `potential :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Output PBE XC potential.
        - `calcEnergy :: LOGICAL, INTENT(IN)`: If `.TRUE.`, calculates the PBE XC energy.
        - `energy :: REAL(KIND=DP), INTENT(OUT)`: Calculated PBE XC energy.
    - **`SUBROUTINE PBE_LibXC(rho, potential, energy)`**:
      - **Description:** Calculates the PBE XC potential and energy using the LibXC library. This routine handles both **spin-unpolarized** (`numSpin=1`) and **spin-polarized** (`numSpin=2`) cases. It initializes PBE exchange (`XC_GGA_X_PBE`) and correlation (`XC_GGA_C_PBE`) functionals from LibXC and then calls `xc_f90_gga_exc_vxc` to get the energy density and potential components. The potential calculation involves terms from `dE/d(rho_up)`, `dE/d(rho_down)`, and `dE/d(sigma_uu)`, `dE/d(sigma_ud)`, `dE/d(sigma_dd)`.
      - **Arguments:**
        - `rho :: REAL(KIND=DP),DIMENSION(:,:,:,:),INTENT(IN)`: Real-space electron density (4th dimension for spin).
        - `potential :: REAL(KIND=DP),DIMENSION(:,:,:,:),INTENT(OUT)`: Output PBE XC potential (spin components).
        - `energy :: REAL(KIND=DP),INTENT(OUT)`: Calculated PBE XC energy.
    - **`FUNCTION PBEStress(rhoReal_SI, rhoRecip_SI) RESULT(stress_tensor)`**:
      - **Description:** Calculates the PBE XC contribution to the stress tensor for **spin-unpolarized** systems using a direct implementation.
      - **Arguments:**
        - `rhoReal_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent real-space density.
        - `rhoRecip_SI :: COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Fourier transform of `rhoReal_SI`.
      - **Return Value:** `stress_tensor :: REAL(KIND=DP), DIMENSION(3,3)`.

# Important Variables/Constants

- **PBE Functional Parameters (defined within routines):**
    - For Exchange: `kappa = 0.8040`, `mu = beta_pbe * (PI^2 / 3.0)` (where `beta_pbe = 0.066724...`).
    - For Correlation: `gamma_pbe = (1.0 - LOG(2.0)) / PI^2`, `beta_pbe` (same as above for correlation gradient term), and parameters `a, a1, b1, b2, b3, b4` for the LDA correlation part within PBE.
- **LibXC Interface (in `PBE_LibXC`):**
    - `x_func, c_func :: TYPE(xc_f90_pointer_t)`: Pointers to LibXC functional objects for PBE exchange and correlation.
    - `XC_GGA_X_PBE, XC_GGA_C_PBE`: Named constants from LibXC identifying the PBE exchange and correlation functionals.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`, `PI`, `IMAG`, `TINY`, `auToGPa`.
- **`PlaneWave`**: For `qVectors` (G-vector components, used for calculating gradients).
- **`FOURIER`**: For the older `FFT` interface (used in `PBEStress` and the direct `PBEPot`).
- **`Fourier_NEW`**: For the newer `FFT_NEW` interface and `FFT_STD_STATE` (used in `PBE_LibXC` and parts of `PBEPot`).
- **`SetupFFT`**: Module is used, likely for global FFT setup consumed by `Fourier_NEW`.
- **`CellInfo`**: For reciprocal space grid dimensions (`k1G, k2G, k3G`), `numSpin`, and `m123G` (total grid points for stress normalization).
- **`MathFunctions`**: `MinMaxVal` is used in `PBEPot`.
- **LibXC Library (External Dependency for `PBE_LibXC`):**
    - `xc_f90_types_m`: For `xc_f90_pointer_t`.
    - `xc_f90_lib_m`: For `xc_f90_func_init`, `xc_f90_gga_exc_vxc`.
    - `libxc_funcs_m`: For functional identifiers like `XC_GGA_X_PBE`.
    This implies that to use `PBE_LibXC`, the code must be compiled and linked against the LibXC library.

**Workflow Selection:**
The choice between the direct PBE implementation (`PBEPot`) and the LibXC version (`PBE_LibXC`) for potential/energy calculation would typically be handled by a higher-level routine (e.g., in `CalPotential`), likely based on `numSpin` (LibXC version is necessary for `numSpin=2`) and possibly compilation flags or input options indicating whether LibXC is available and preferred. The `exchangeCorrelation` flag (from `XC_LDA` module, though its name is generic) would select PBE (value 2) over LDA.
