# Overview

The `KEDF_HC10` module implements the nonlocal orbital-free Kinetic Energy Density Functional (KEDF) developed by Huang and Carter, specifically tailored for semiconductors, as described in Phys. Rev. B 81, 045206 (2010). This functional aims to improve KEDF accuracy by using a density-dependent kernel.

The core of the HC10 KEDF involves an energy expression that depends on integrals of `rho^alpha(r) * w(r-r') * rho^beta(r')`, where `w` is a kernel function. The kernel itself is not fixed but depends on a local variable `xi(r) = kF(r) * (1 + lambda_hc * s(r)^2)`, where `kF` is the local Fermi wavevector and `s` is the reduced density gradient. The actual calculation involves interpolating kernel values (or their Fourier transforms) based on a set of reference `xi` points.

This module provides functions to calculate the HC10 kinetic energy (`intEnergy`), its corresponding potential (`intPot`), and its contribution to the stress tensor (`intStress`).

# Key Components

- **`MODULE KEDF_HC10`**: The main container for HC10 KEDF functionalities.

- **`SUBROUTINE MakeHC10Kernel(xi_val, k0, k1)`**:
  - **Description:** Computes the Fourier transforms of the fundamental HC10 kernel components for a given reference `xi_val`. Specifically, `k0` becomes the Fourier transform of `w(eta(G; xi_val))` and `k1` becomes the Fourier transform of `eta * dw/deta`. The underlying universal kernel `w(eta)` and its derivative `wp(eta)` (where `eta = G / (2*xi_val)`) are obtained by solving an ODE via the `IntKernelODE` module (specifically `makeKernel_ode1`), if not already done.
  - **Arguments:**
    - `xi_val :: REAL(KIND=DP), INTENT(IN)`: The reference `xi` value for which the kernel transforms are computed.
    - `k0, k1 :: COMPLEX(KIND=DP), INTENT(OUT), DIMENSION(:,:,:)`: Output arrays for the Fourier transformed kernel components.

- **`FUNCTION intEnergy(rho) RESULT(energy_val)`**:
  - **Description:** Calculates the HC10 kinetic energy. It first computes `xi(r)` for the input density `rho(r)`. Then, it iterates through bins defined by `refxi` values. For each bin, it uses cubic Hermite interpolation (functions `h00`, `h01`, `h10`, `h11`) with pre-computed kernel transforms (`kernel0`, `kernel1` from `MakeHC10Kernel` at bin boundaries) to construct the effective G-space kernel for that density region. The energy is then computed as `cTF * C * SUM(rho**beta * FFT(FFTx))`, where `FFTx` is the sum of interpolated kernel contributions multiplied by `FFT(rhoMask*h_interp_func)`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION intPot(rho, energy) RESULT(potential_val)`**:
  - **Description:** Calculates the potential corresponding to `intEnergy`. The calculation involves similar Hermite interpolation steps for different terms arising from the functional derivative of the HC10 energy expression. It computes contributions from `d(rho^alpha)/drho`, `d(rho^beta)/drho`, and derivatives of `xi` with respect to `rho` and `grad(rho)`.
  - **Arguments:**
    - `rho :: REAL(KIND=DP), DIMENSION(n1G,n2G,n3G), INTENT(IN)`.
    - `energy :: REAL(KIND=DP), INTENT(OUT)`: The calculated HC10 energy (recalculated within this function).

- **`FUNCTION intStress(vol, rho) RESULT(stress_tensor)`**:
  - **Description:** Calculates the HC10 KEDF contribution to the stress tensor. This is a complex calculation involving derivatives of the kernel with respect to `eta` and further interpolation steps, similar to `intEnergy` and `intPot`.
  - **Arguments:**
    - `vol :: REAL(KIND=DP), INTENT(IN)`: Cell volume.
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`SUBROUTINE MakeXi(rho, xi, ss)`**:
  - **Description:** Calculates the local variable `xi(r) = kF(r) * (1 + hc_lambda_val * ss(r)^2)` and the reduced gradient squared `ss(r) = |nabla rho(r)|^2 / rho(r)^(8/3)`. This version calculates `nabla rho` internally using FFTs.
  - **Arguments:** `rho`, `xi`, `ss` (outputs).

- **`SUBROUTINE MakeXiG(rho, xi, ss, gradRho, rho83)`**:
  - **Description:** Similar to `MakeXi`, but takes pre-calculated `gradRho` and `rho^(8/3)` (`rho83`) as input, avoiding redundant FFTs.
  - **Arguments:** `rho`, `gradRho`, `rho83` (inputs); `xi`, `ss` (outputs).

- **`FUNCTION h00(t)`, `h10(t)`, `h01(t)`, `h11(t)`**:
  - **Description:** Element-wise functions that return the values of the cubic Hermite spline basis functions for a given parameter `t` (where `t` is typically `(xi(r) - xi_ref_lower) / (xi_ref_upper - xi_ref_lower)`).

# Important Variables/Constants

- **Module-Level Parameters for HC10 KEDF:**
    - `hc_lambda_val :: REAL(KIND=DP)`: The `lambda` parameter in the definition of `xi` (default: 0.01).
    - `cutrho :: REAL(KIND=DP)`: A density cutoff below which contributions might be ignored (default: 0.0).
    - `refRatio :: REAL(KIND=DP)`: Ratio used to generate a geometric series of reference `xi` points (`refxi`) for interpolation (default: 1.15).
    - `refxi(:) :: REAL(KIND=DP), ALLOCATABLE`: Array of reference `xi` values.
    - `nref :: INTEGER`: Number of points in `refxi`, determined dynamically based on `min(xi)` and `max(xi)` in the system.
    - `have_RK :: LOGICAL`: Flag (default `.FALSE.`) set to `.TRUE.` after the fundamental kernel `w(eta)` is first computed by `IntKernelODE` via `MakeHC10Kernel`.
- **From `KEDF_WTkernel` (via `USE` statement):**
    - `alpha, beta :: REAL(KIND=DP)`: Exponents used for `rho^alpha` and `rho^beta` terms in the energy expression.
- **From `KEDF_TF` (via `USE` statement):**
    - `cTF :: REAL(KIND=DP)`: The Thomas-Fermi constant `(3/10)*(3*pi^2)^(2/3)`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`, `PI`, `IMAG`.
- **`CellInfo`**: For grid dimensions (`k1G`, `k2G`, `k3G`, `k123G`, `n1G`, `n2G`, `n3G`, `n123G`) and `cell` variable (for `cell%dV` in stress).
- **`Sys`**: For `rhoS` (reference density used for generating `refxi` scale).
- **`KEDF_WTkernel`**: For `alpha` and `beta` parameters. This suggests HC10 might be related to or share formalism with the Wang-Teter KEDF.
- **`PlaneWave`**: For `qTable` (G-vector magnitudes), `qMask`, `qVectors`.
- **`FOURIER`**: For the `FFT` interface, used extensively for convolutions and gradient calculations.
- **`MPI_Functions`**: For parallel reductions (`MPI_ALLREDUCE`) and utilities (`TITLE`, `StartClock`, `StopClock`).
- **`OUTPUT`**: For `WrtOut`.
- **`KEDF_TF`**: For `cTF`.
- **`MathFunctions`**: For `GPrime` (derivative of Lindhard function for stress) and `stepfun` (step function for masking).
- **`IntKernelODE`**: This is a critical dependency for `MakeHC10Kernel`. It provides:
    - `ode_beta` (set from HC10's `beta`).
    - `wInf` (asymptotic value for kernel ODE).
    - `int_eta`, `int_w`, `int_wp` (results from ODE solution: `eta` grid, `w(eta)`, and `eta*w'(eta)`).
    - `makeKernel_ode1` (the ODE solver).
    - `IntKernelODEClean` (aliased as `Clean`).

The HC10 functional is used by first calling `MakeHC10Kernel` (which itself calls `IntKernelODE` routines if `have_RK` is false) to tabulate the base kernel `w(eta)`. Then, functions like `intEnergy` and `intPot` use these tabulated values, along with Hermite interpolation based on the local `xi` variable, to calculate the energy and potential. The `xi` variable introduces a dependency on both `rho` and `nabla rho`.
