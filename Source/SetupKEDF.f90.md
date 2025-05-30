# Overview

The `SetupKEDF` module plays a crucial role in initializing the parameters and data structures required for the chosen Kinetic Energy Density Functional (KEDF) in a PROFESS simulation. Based on an integer flag (`kinetic` from the `CellInfo` module), it configures various KEDF modules by setting their specific parameters (like `lambda` for TF, `mu` for vW, `alpha`, `beta`, `gamma` for nonlocal functionals, etc.) and reference densities (`rho0`, `rhoS`).

For nonlocal KEDFs (WT, WGC, HC10, CAT, MGP), this module is responsible for allocating the shared `keKernel` array (which resides in `KEDF_WTkernel`) and then calling the appropriate routine (e.g., `FillWT`, `FillWGC`, `FillCAT`, `FillMGP`) to populate this kernel with values specific to the chosen functional.

Additionally, if Particle Mesh Ewald (PME) methods are enabled for ion-electron or ion-ion interactions (via `ieSpline` or `iiSpline` flags), `SetupFunctional` also handles the setup of B-spline coefficients by calling `FillB` from the `CBSpline` module.

The `KEDFRefresh` subroutine is designed to be called whenever the simulation cell changes, to update G-space dependent quantities like `qTable` and to regenerate the nonlocal kernels if necessary.

# Key Components

- **`MODULE SetupKEDF`**: The main container module.

- **`SUBROUTINE SetupFunctional`**:
  - **Description:** This is the primary initialization routine for all KEDF-related settings.
    1.  Allocates G-space arrays: `qTable`, `qVectors`, `qMask` (from `PlaneWave` module).
    2.  If PME is active (`ieSpline` or `iiSpline` from `IonElectronSpline` is true):
        - Allocates `bSpline1`, `bSpline2`, `bSpline3` (from `CBSpline`).
        - Sets/checks `splineOrder`.
        - Calls `FillB` to compute B-spline coefficients.
    3.  Uses a `SELECT CASE(kinetic)` block to perform specific setups for the chosen KEDF:
        - **Case 1 (TF):** Sets `lambda` (in `KEDF_TF`), `mu=0`.
        - **Case 2 (VW):** Sets `mu` (in `KEDF_VW`), `lambda=0`.
        - **Case 3 (TF+VW):** Sets `lambda` and `mu`.
        - **Case 4 (WT), 17 (EvW-WT):** Sets `lambda`, `mu`, `alpha`, `beta` (in `KEDF_WTkernel`), `rho0` (in `SYS`). Allocates `keKernel(k1G,k2G,k3G,1)`. Sets `kloc`, `aloc` (for EvW, in `KEDF_EvW`).
        - **Case 5 (WGC), 12 (WGCD), 18 (EvW-WGC):** Sets `lambda`, `mu`. Sets WGC parameters `alpha`, `beta`, `gamma`, `alpha5`, `beta5` (in `KEDF_WGCkernel` and `KEDF_WTkernel`). Sets reference densities `rho0`, `rhoS` and flags `hold0`, `holdS` (in `SYS`). Allocates `keKernel(k1G,k2G,k3G,4)`. For WGCD (case 12), also initializes `t`, `scalefunt`, `scalefunDD` arrays (in `KEDF_WGCD`) for spline interpolation of the F(r) function. Sets `kloc`, `aloc` (for EvW).
        - **Case 7 (LQ), 8 (HQ):** Sets `lambda=1`, `mu=1`. Sets `rho0`.
        - **Case 10 (CAT):** Allocates `keKernel(k1G,k2G,k3G,4)`. Sets `rho0`, `rhoS`.
        - **Case 11 (HC10):** Allocates `hc_lambda` array (in `KEDF_HC10`). Sets `lambda=1`, `mu=1`. Sets `hc_lambda_val`, `alpha`, `beta`, `rhoS`, `rho0`.
        - **Case 15 (GGA):** Sets `lambda`, `mu`. Initializes GGA parameters in `KEDF_GGA` module (e.g., `GGA_functional`, `Lmu`, etc.).
        - **Case 16 (GGA+WGCD):** Sets up `rho0` and interpolation table for ELF in `KEDF_GGA`.
        - **Case 19 (MGP):** Allocates `keKernel(k1G,k2G,k3G,3)`. Sets `lambda`, `mu`, `rho0`, `gamma`. Initializes `rhoS`.
        - **Default:** Error and stop.
  - **Side Effects:** Modifies parameters in many external KEDF modules, allocates several key arrays.

- **`SUBROUTINE KEDFRefresh(kinetic)`**:
  - **Description:** This routine is called when the simulation cell changes. It first calls `FillQTable` (from `PlaneWave`) to update the G-vector table (`qTable`, `qVectors`, `qMask`) based on the new `cell%cellRecip`. Then, based on the `kinetic` flag, it calls the appropriate kernel-filling routine for nonlocal KEDFs (`FillWT`, `FillWGC`, `FillMGP`, `FillCAT`), after potentially updating reference densities `rho0` and `rhoS` if they are not held constant (`hold0`, `holdS` flags).
  - **Arguments:** `kinetic :: INTEGER, INTENT(IN)`.

# Important Variables/Constants

This module primarily *sets* variables in other modules. Key input parameters that control its behavior:
- `kinetic :: INTEGER` (from `CellInfo`): Determines which KEDF pathway is configured.
- `energyCutoff :: REAL(KIND=DP)` (passed to `InitializeFFT` which is called by `Sys::SetupSystem` before `SetupFunctional`): Implicitly determines FFT grid sizes (`k1G, k2G, k3G`) which are then used for allocating kernels.
- Various parameters from specific KEDF modules are initialized to default values if they haven't been set (e.g. `lambda`, `mu`, `alpha`, `beta`, `gamma`, `rho0`, `rhoS`).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

`SetupKEDF` is a high-level orchestrator with extensive dependencies:
- **Core Utilities:** `CONSTANTS`, `MPI_Functions`, `Output`, `OutputFiles`.
- **System & Cell Data:** `SYS` (for setting `rho0`, `rhoS`, `hold0`, `holdS`), `CellInfo` (for `kinetic` flag, `cell` data, grid dimensions, `numSpin`).
- **FFT & Plane Waves:** `PlaneWave` (for allocating and filling `qTable`, `qVectors`, `qMask` via `FillQTable` in `KEDFRefresh`).
- **Splines (for PME):** `IonElectronSpline` (for `iiSpline`, `ieSpline` flags), `CBSpline` (for `FillB`, `splineOrder`, and `bSpline` arrays).
- **Specific KEDF Modules (for setting parameters and calling kernel fillers):**
    - `KEDF_TF` (sets `lambda`)
    - `KEDF_VW` (sets `mu`)
    - `KEDF_WTkernel` (sets `alpha`, `beta`; allocates `keKernel`; `FillWT` is called by `KEDFRefresh`)
    - `KEDF_WGCkernel` (sets `firstOrderWGC`, `gamma`, `alpha5`, `beta5`; `FillWGC` is called by `KEDFRefresh`)
    - `KEDF_HC10` (sets `hc_lambda_val` and allocates `hc_lambda`)
    - `KEDF_WGCD` (initializes `t`, `scalefunt`, `scalefunDD` arrays; uses `mrhos`)
    - `KEDF_EvW` (sets `kloc`, `aloc`)
    - `KEDF_GGA` (initializes various parameters for GGA functionals)
    - `KEDF_MGPkernel` (`FillMGP` is called by `KEDFRefresh`)
    - `KEDF_CAT` (`FillCAT` is called by `KEDFRefresh`)
- **Math Utilities:** `MathSplines` (used by `KEDF_WGCD` for its scaling function setup).

The `SetupFunctional` routine is typically called once during the main program initialization (e.g., from `Initializer` module). `KEDFRefresh` is called whenever the cell geometry changes, for example, during cell optimization.
