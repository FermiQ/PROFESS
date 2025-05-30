# Overview

The `KEDF_EvW` module implements the Embedded von Weizsacker (EvW) Kinetic Energy Density Functional (KEDF). This KEDF is particularly designed for semiconductors and aims to improve KEDF accuracy by enhancing the von Weizsacker (vW) term. The enhancement factor, `aloc`, is related to the localization of the electron density and is determined self-consistently in a more complete implementation (though in the provided routines, `aloc` is used as a fixed parameter).

The EvW-KEDF typically modifies existing KEDFs (like Thomas-Fermi + Wang-Teter or Thomas-Fermi + Wang-Govind-Carter) by scaling their components:
- The Thomas-Fermi (TF) and the non-local (WT or WGC) parts are scaled by `(1 - aloc/2)`.
- The von Weizsacker (vW) part is scaled by `(1 + aloc)`.

This module provides two main subroutines: `Cal_EVT` (EvW based on Wang-Teter) and `Cal_EVC` (EvW based on Wang-Govind-Carter).

**Reference:** J. Chem. Phys. 140, 18A531 (2014) Shin and E. A. Carter.

# Key Components

- **`MODULE KEDF_EvW`**: The main container for the EvW-KEDF calculations.

- **`SUBROUTINE Cal_EVT(potential, rho, calcEnergy, locETable, optSqrt)`**:
  - **Description:** Calculates the EvW-KEDF potential and energy, where the non-local component is based on the Wang-Teter (WT) KEDF.
    The total KEDF energy is `(1-aloc/2) * (E_TF + E_WT) + (1+aloc) * E_vW`.
    The potential is derived accordingly.
  - **Arguments:**
    - `potential :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT)`: Output KEDF potential.
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: Input electron density (assumed spin-unpolarized for EvW, uses `rho(:,:,:,1)`).
    - `calcEnergy :: LOGICAL, INTENT(IN)`: Flag to compute energy components.
    - `locETable :: REAL(KIND=DP), DIMENSION(9), INTENT(OUT)`: Array to store energy components (`locETable(7)` for TF, `locETable(8)` for vW, `locETable(9)` for WT).
    - `optSqrt :: LOGICAL, INTENT(IN)`: Flag indicating if optimization is w.r.t. `sqrt(rho)`.

- **`SUBROUTINE Cal_EVC(potential, rho, calcEnergy, locETable, optSqrt)`**:
  - **Description:** Calculates the EvW-KEDF potential and energy, where the non-local component is based on the Wang-Govind-Carter (WGC) KEDF.
    The total KEDF energy is `(1-aloc/2) * (E_TF + E_WGC) + (1+aloc) * E_vW`.
    The potential is derived accordingly.
  - **Arguments:** Same as `Cal_EVT`, but `locETable(9)` will store the WGC energy component.

# Important Variables/Constants

- **Module-Level Parameters for EvW-KEDF:**
    - **`kloc :: REAL(KIND=DP)`**: A parameter for the EvW KEDF (default: -100.0). Its direct usage is not apparent in `Cal_EVT` or `Cal_EVC`, suggesting it might be used in the self-consistent determination of `aloc` or other parts of the EvW scheme not shown.
    - **`aloc :: REAL(KIND=DP)`**: The primary parameter for the EvW KEDF (default: -100.0). It controls the scaling of the TF/non-local and vW components of the kinetic energy and potential. In a full EvW implementation, `aloc` would be determined self-consistently based on density localization. In these routines, it acts as a fixed input parameter.
    - **`tolk :: REAL(KIND=DP)`**: Tolerance for the iterative self-consistent loop that would typically determine `aloc` (default: 1.0e-8). This parameter is not used by the provided subroutines as they don't perform the `aloc` convergence loop.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP` (double precision kind parameter).
- **`MPI_Functions`**: For `Error` and `message` (error handling and messaging). `TITLE` is imported but not used in these subroutines.
- **`KEDF_TF`**: For `TFPotential` and `TFEnergy` (calculates Thomas-Fermi contributions).
- **`KEDF_VW`**: For `VWPotentialSqrtPlus` (calculates von Weizsacker contributions, taking `sqrt(rho)` as input).
- **`CellInfo`**: For `cell` (cell information), `numSpin` (number of spin channels; EvW routines currently stop if `numSpin == 2`), and grid dimensions (`n1G`, `n2G`, `n3G`).
- **`SYS`**: For `bvac` (boolean flag indicating vacuum treatment, passed to WT/WGC potential routines).
- **`KEDF_WT`**: The `Cal_EVT` subroutine uses `WTPotentialPlus` to get the Wang-Teter potential and energy.
- **`KEDF_WGC`**: The `Cal_EVC` subroutine uses `WGCPotentialPlus` to get the Wang-Govind-Carter potential and energy.

**Workflow:**
When `kinetic` flag in `CalPotPlus` selects EvW-KEDF (e.g., type 17 for EvW-WT, type 18 for EvW-WGC), the corresponding subroutine (`Cal_EVT` or `Cal_EVC`) is called. These routines:
1. Check if the calculation is spin-polarized (currently, EvW is implemented for `numSpin=1`).
2. Calculate the TF contribution, scaled by `(1-aloc/2)`.
3. Calculate the non-local (WT or WGC) contribution, also scaled by `(1-aloc/2)`.
4. If optimizing w.r.t `sqrt(rho)`, apply the chain rule to the sum of TF and non-local potentials.
5. Calculate the vW contribution, scaled by `(1+aloc)`.
6. Combine these components to form the total EvW KEDF potential.
7. If `calcEnergy` is true, the corresponding scaled energy terms are stored in `locETable`.

The determination of the `aloc` parameter itself is not handled within these specific subroutines; they assume `aloc` is a pre-set module-level variable. A complete EvW implementation would involve an outer loop to determine `aloc` self-consistently.
