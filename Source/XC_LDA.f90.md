# Overview

The `XC_LDA` module implements the Local Density Approximation (LDA) and Local Spin Density Approximation (LSDA) for the exchange-correlation (XC) energy and potential in Density Functional Theory (DFT). It provides functions to calculate these quantities as well as their contributions to the stress tensor.

For the spin-unpolarized case (LDA), the Perdew-Zunger (PZ) parametrization (Phys. Rev. B 23, 5048 (1981)) is primarily used. For the spin-polarized case (LSDA), the module offers:
1.  Perdew-Zunger (PZ) LSDA.
2.  Perdew-Wang 1992 (PW92) LSDA, which uses the `spn_gcor` helper subroutine for its correlation part.

The choice between these is typically controlled by the module-level `lsda` parameter. The overall selection of LDA/LSDA versus other functionals (like GGA) is controlled by the `exchangeCorrelation` parameter.

# Key Components

- **`MODULE XC_LDA`**: The main container for LDA/LSDA XC calculations.

- **Module-Level Parameters:**
    - `exchangeCorrelation :: INTEGER`: A global flag (often set from input) indicating the type of XC functional to use. A value of `1` typically selects LDA/LSDA, meaning routines from this module will be called. (Default: 1)
    - `lsda :: INTEGER`: A flag to select the specific LSDA parametrization when spin-polarized calculations are performed (`numSpin = 2`).
        - `1`: Perdew-Zunger (PZ) LSDA. (Default)
        - `2`: Perdew-Wang 1992 (PW92) LSDA.

- **Energy and Potential Functions/Subroutines:**
    - **`FUNCTION LDAEnergy(rhoReal) RESULT(energy_val)`**:
      - **Description:** Calculates the XC energy.
        - If `SIZE(rhoReal,4) == 1` (spin-unpolarized): Uses the Perdew-Zunger LDA formula for exchange and correlation.
        - If `SIZE(rhoReal,4) == 2` (spin-polarized): Uses the Perdew-Wang 92 LSDA formula (via `spn_gcor` for correlation) for exchange and correlation.
      - **Arguments:** `rhoReal :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)` (Real-space density, 4th dim is spin).
    - **`FUNCTION LDAPot(rho) RESULT(potential_val)`**:
      - **Description:** Calculates the LDA XC potential for spin-unpolarized systems using the Perdew-Zunger parametrization.
      - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)` (expects 4th dim to be 1).
      - **Return Value:** `potential_val :: REAL(KIND=DP), DIMENSION(SIZE(rho)...)`
    - **`SUBROUTINE LSDAPotPW92(rho, LDAPotential, LDAEnergy)`**:
      - **Description:** Calculates the LSDA XC potential and energy using the Perdew-Wang 92 parametrization for spin-polarized systems. It calls the `spn_gcor` helper for the correlation part.
      - **Arguments:**
        - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)` (4th dim is 2).
        - `LDAPotential :: REAL(KIND=DP), DIMENSION(SIZE(rho)...), INTENT(OUT)`.
        - `LDAEnergy :: REAL(KIND=DP), INTENT(OUT)`.
    - **`SUBROUTINE LSDAPotPZ(rho, LDAPotential, LDAEnergy)`**:
      - **Description:** Calculates the LSDA XC potential and energy using the Perdew-Zunger parametrization for spin-polarized systems.
      - **Arguments:** Same as `LSDAPotPW92`.
    - **`SUBROUTINE spn_gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)`**:
      - **Description:** A helper subroutine that computes components of the local spin correlation energy (`gg`) and its derivative with respect to `rs` (`ggrs`) based on the Perdew-Wang 92 parametrization for given parameters `a, a1, ...` and `rs`.
      - **Arguments:** Various PW92 parameters, `rs` (IN), `gg` (OUT), `ggrs` (OUT).

- **Stress Functions:**
    - **`FUNCTION LSDAStress(cellVol, rhoR, energy) RESULT(stress_tensor)`**:
      - **Description:** Calculates the LSDA contribution to the stress tensor for spin-polarized systems. It uses `LSDAPotPZ` internally to get the potential needed for the stress formula `Stress_aa = E_xc/Vol - SUM_spin SUM_grid { rho_spin * V_xc_spin } / N_grid_total`.
      - **Arguments:** `cellVol`, `rhoR` (4D real-space density), `energy` (LSDA XC energy).
    - **`FUNCTION LDAStress(cellVol, rhoR, energy) RESULT(stress_tensor)`**:
      - **Description:** Calculates the LDA contribution to the stress tensor for spin-unpolarized systems. Uses `LDAPot`.
      - **Arguments:** `cellVol`, `rhoR` (4D real-space density, 4th dim is 1), `energy` (LDA XC energy).

# Important Variables/Constants

- **Perdew-Zunger (PZ) Parameters:**
    - `cX = -0.73855876638202234_DP`: Exchange constant `(-3/4) * (3/pi)^(1/3)`.
    - `cC = 0.62035049089940009_DP`: `(3/(4*pi))^(1/3)`.
    - `a(2), b(2), c(2), d(2)`: Parameters for PZ correlation at `rs < 1` (paramagnetic and ferromagnetic).
    - `g(2), b1(2), b2(2)`: Parameters for PZ correlation at `rs >= 1`.
    - Constants `cC1-cC5` in `LDAPot` are derived from these for potential calculation.
- **Perdew-Wang 92 (PW92) Parameters (used in `spn_gcor` via `LSDAPotPW92` and `LDAEnergy`):**
    - Parameters like `0.0310907d0, 0.21370d0, ...` are passed to `spn_gcor`.
    - `fzz = 1.709921d0`, `gamma = 0.5198421d0` are spin-scaling interpolation parameters.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision) and `PI`.
- **`CellInfo`**: The stress functions (`LDAStress`, `LSDAStress`) import `m123G` (total global grid points), `m3G` (global Z dimension), and `n3G` (local Z dimension for current processor) for correct normalization and scaling of stress contributions in parallel execution.
- **`Timer`**: (Implicit) Calls to `StartClock` and `StopClock` are present in `LSDAPotPW92` and `LSDAPotPZ`, suggesting use of the `Timer` module for profiling.

**Workflow:**
Typically, a higher-level routine (e.g., in `CalPotential` module) would:
1. Check the global `exchangeCorrelation` flag. If it's 1 (for LDA/LSDA):
2. Check `numSpin` (from `CellInfo`).
3. If `numSpin == 1`, call `LDAEnergy` and `LDAPot`.
4. If `numSpin == 2`, check the `lsda` flag:
   - If `lsda == 1` (PZ), call `LSDAPotPZ` (which provides both energy and potential).
   - If `lsda == 2` (PW92), call `LSDAPotPW92` (which provides both energy and potential).
The stress functions are called similarly by routines in `CalStress`.
