# Overview

The `KEDF_GGA` module implements a variety of Generalized Gradient Approximation (GGA) and some Laplacian-level meta-GGA Kinetic Energy Density Functionals (KEDFs). GGA KEDFs improve upon the local Thomas-Fermi (TF) functional by incorporating information about the gradient of the electron density, typically through a dimensionless reduced gradient `s`. Laplacian-level meta-GGAs further include the Laplacian of the density, often via a reduced Laplacian `q`.

This module provides two main subroutines:
1.  `GGAPotentialPlus`: A versatile routine that calculates the kinetic energy and potential for a wide range of predefined GGA and Laplacian-level meta-GGA KEDFs. The specific functional is selected by an integer flag.
2.  `vWGTF`: Implements a specific KEDF model of the form `T_vW + G(rho) * T_TF`, where `T_vW` is the von Weizsacker energy, `T_TF` is the Thomas-Fermi energy, and `G(rho)` is a density-dependent scaling factor. Different forms of `G(rho)` can be chosen.

**References for `GGAPotentialPlus`:**
- JCTC 5, 3161 (2009)
- JCP 127, 144109 (2007)
- Laplacian-Level Meta-GGA: L.A. Constantin, E. Fabiano and F. Della Sala

**Reference for `vWGTF`:**
- J. Xia and E. A. Carter, submitted to PRB (at the time of comment)

# Key Components

- **`MODULE KEDF_GGA`**: The main container module.

- **`SUBROUTINE GGAPotentialPlus(rho, potential, calcEnergy, energy, functionals)`**:
  - **Description:** Calculates the kinetic energy and potential for various GGA and Laplacian-level meta-GGA KEDFs. The core of these functionals is an enhancement factor `F(s)` (for GGAs) or `F(s,q)` (for meta-GGAs) which multiplies the Thomas-Fermi energy density. The subroutine computes `s = |nabla rho| / (cs * rho^(4/3))` and `q = laplacian(rho) / (cs^2 * rho^(5/3))`, then evaluates `F` and its derivatives `dF/ds`, `dF/dq` based on the integer `functionals` flag.
  - **Arguments:**
    - `rho :: REAL(kind=DP), Dimension(:,:,:), INTENT(in)`: Total electron density.
    - `potential :: REAL(kind=DP), Dimension(:,:,:), INTENT(out)`: Output KEDF potential (`dE_k/drho`).
    - `calcEnergy :: LOGICAL, INTENT(in)`: If true, calculates the kinetic energy.
    - `energy :: REAL(kind=DP), INTENT(out)`: Calculated kinetic energy.
    - `functionals :: INTEGER, INTENT(in)`: Flag selecting the specific GGA/meta-GGA functional.
      - `1`: TF only (`F=1`)
      - `2`: vW only (`F = s^2/(8*cTF)`)
      - `3`: `lambda*TF + mu*vW`
      - Other cases implement various published GGA/meta-GGA forms (LLP, OL2, PW91k, TW02, PBE2, LC94, E00, P92, PW86, DK, OL1, B86A, B86B, Thak, DK87, LKT, GE4, PGSL, PGSLr).
      - `functionals >= 50` are Laplacian-level meta-GGAs.

- **`SUBROUTINE vWGTF(rho, potential, calcEnergy, energy, model)`**:
  - **Description:** Calculates kinetic energy and potential for a model of the form `T[rho] = T_vW[rho] + G(rho) * T_TF[rho]`. `T_vW` is the von Weizsacker term. `G(rho)` is a scaling factor that depends on the density `rho` relative to a reference density `rho0`.
  - **Arguments:**
    - `rho :: REAL(kind=DP), Dimension(:,:,:), INTENT(in)`: Total electron density.
    - `potential :: REAL(kind=DP), Dimension(:,:,:), INTENT(out)`: Output KEDF potential.
    - `calcEnergy :: LOGICAL, INTENT(in)`: If true, calculates the kinetic energy.
    - `energy :: REAL(kind=DP), INTENT(out)`: Calculated kinetic energy.
    - `model :: INTEGER, INTENT(in)`: Selects the form of `G(rho)`:
      - `1`: `G = co * (rho/rho0)^expo`
      - `2`: `G = sqrt(1/ELF - 1)`, where `ELF` is a function of `(rho/rho0)^gb`.
      - `3`: `G` derived from an interpolated Electron Localization Function (ELF).

# Important Variables/Constants

- **Module-Level Parameters:**
    - `GGA_functional :: INTEGER`: Default functional selector for `GGAPotentialPlus` (default: 1, TF only).
    - `model :: INTEGER`: Default model selector for `vWGTF` (default: 1).
    - `CP :: REAL(DP)`: Penalty coefficient for an optional vW term in `GGAPotentialPlus` (default: 0.0). If non-zero, a penalty energy related to the difference between exact vW and its gradient expansion is added.
    - `LKTa0, Lmu, Lbet, Llam, Lsig`: Parameters for specific GGA/meta-GGA functionals (LKT, PGSLr).
    - `numint, d, ELFvsd, ELFvsdDD`: Parameters and arrays for spline interpolation in `vWGTF` (model 3).

- **Key Derived Quantities in `GGAPotentialPlus`:**
    - `s :: REAL(KIND=DP), Dimension(:,:,:)`: Reduced density gradient.
    - `qq :: REAL(KIND=DP), Dimension(:,:,:)`: Reduced Laplacian of the density (for `functionals >= 50`).
    - `F :: REAL(KIND=DP), Dimension(:,:,:)`: The GGA/meta-GGA enhancement factor.
    - `dFds, dFdqq`: Derivatives of `F` w.r.t. `s` and `qq`.

- **Key Derived Quantities in `vWGTF`:**
    - `G :: REAL(KIND=DP), Dimension(:,:,:)`: The scaling factor for the TF term.
    - `dGdrho :: REAL(KIND=DP), Dimension(:,:,:)`: Derivative of `G` w.r.t. `rho`.

- **Constants for Functionals:** Numerous `PARAMETER` constants (e.g., `cTF`, `cs`, `LLPb`, `PWA1`, etc.) define the specific mathematical forms of the various GGA functionals implemented in `GGAPotentialPlus`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP` (double precision), `pi`, `imag`.
- **`Fourier`**: For the `FFT` interface, used to calculate density gradients and Laplacians via reciprocal space.
- **`OutputFiles`**: For `outputUnit` (standard output).
- **`KEDF_TF`**: For `lambda` (TF coefficient, imported but not directly used by these routines, TF is handled internally).
- **`KEDF_VW`**: For `mu` (vW coefficient, imported but not directly used, vW handled internally or in `vWGTF`).
- **`PlaneWave`**: For `qTable` (G^2 values), `qVectors` (G-vector components), `qMask` (G-vector mask). These are used when calculating gradients/Laplacians.
- **`CellInfo`**: For reciprocal space grid dimensions (`k1G`, `k2G`, `k3G`).
- **`MathSplines`**: Used by `vWGTF` (model 3) for `spline_cubic_set` and `spline_cubic_val`.
- **`Sys`**: Used by `vWGTF` for `rho0` (reference density).

The choice of which GGA functional or `vWGTF` model to use is typically controlled by input parameters that set the `GGA_functional` or `model` integer flags at the module level or are passed into the routines. The `GGAPotentialPlus` routine is a large `SELECT CASE` block that implements the mathematical formulas for many different KEDFs.
