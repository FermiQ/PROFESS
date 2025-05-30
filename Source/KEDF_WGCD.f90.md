# Overview

The `KEDF_WGCD` module implements the Wang-Govind-Carter with Density Decomposition (WGCD) Kinetic Energy Density Functional (KEDF). This advanced KEDF is based on the work by Junchao Xia and Emily A. Carter, "Density-decomposed orbital-free density functional theory for covalently-bonded molecules and materials," Phys. Rev. B 86, 235109 (2012).

The core idea of WGCD is to decompose the total electron density `rho(r)` into a smooth part `rho1(r)` and a localized part `rho2(r)`. This decomposition is achieved using a spatially varying scaling function `F(r)` such that `rho1(r) = rho(r) * F(r)` and `rho2(r) = rho(r) * (1 - F(r))`. The scaling function `F(r)` itself depends on the local density relative to a reference density `rho0` (i.e., `F(r) = f(rho(r)/rho0)`), and is determined self-consistently.

The kinetic energy is then formulated as:
`T_s[rho] = (T_1[rho] - T_1[rho1]) + T_2[rho1]`
where:
- `T_1[rho]` is typically a simpler KEDF, often `lambda*TF + mu*vW` (Thomas-Fermi + von Weizsacker).
- `T_2[rho1]` is a more sophisticated KEDF applied to the smooth density `rho1`, usually `TF[rho1] + vW[rho1] + WGC[rho1]` (Wang-Govind-Carter).

The module provides routines to initialize this iterative process, compute the `F(r)` matrix self-consistently, and then calculate the total WGCD kinetic energy and potential.

# Key Components

- **`MODULE KEDF_WGCD`**: The main container for WGCD KEDF functionalities.

- **`SUBROUTINE InitializeWGCD(dimX, dimY, dimZ, kinetic, pot_tol)`**:
  - **Description:** Initializes the WGCD calculation. It allocates space for `Fr` (the scaling function `F(r)`) and `oldFr`. Crucially, it sets `Fr = 1.0` for the first iteration, meaning `rho1 = rho` initially. It temporarily modifies the `kinetic` KEDF type to 4 (TF+vW+WT, Wang-Teter) and adjusts `alpha` and `beta` parameters (from `KEDF_WTkernel`) to ensure a stable start for determining an initial `rho`. The `pot_tol` for the electron density optimization is also temporarily relaxed.
  - **Arguments:**
    - `dimX, dimY, dimZ :: INTEGER, INTENT(IN)`: Dimensions for allocating `Fr`.
    - `kinetic :: INTEGER, INTENT(INOUT)`: KEDF type selector, modified to 4.
    - `pot_tol :: REAL(KIND=DP), INTENT(INOUT)`: Potential tolerance for SCF, temporarily increased.

- **`SUBROUTINE ComputeFrMatrix(rho, dimX, dimY, dimZ, nspin, new_iteration, kinetic, pot_tol)`**:
  - **Description:** This is the core of the self-consistent determination of `Fr`.
    1.  On the very first call after `InitializeWGCD` (when `WGCDiter==0`), it switches `kinetic` to 12 (the WGCD KEDF type, which internally will use WGC for `T_2`), restores `alpha` and `beta` to their WGC values, and reallocates `keKernel` for WGC.
    2.  It increments `WGCDiter`.
    3.  `oldFr` is updated with the current `Fr`.
    4.  A new `rho0` is calculated based on `SUM(rho*Fr) / SUM(Fr_mask_for_rho0)`, where `Fr_mask_for_rho0` might exclude vacuum regions if `bvac` is true.
    5.  The new `Fr(r)` is computed by interpolating a pre-defined scaling function `scalefunt` using `rho(r)/rho0 - shiftm` as the argument. The interpolation uses `spline_cubic_val`.
    6.  A new `rhoS` (reference density for WGC kernel) is determined from the new `rho0`.
    7.  The WGC kernel is refilled using `FillWGC` with the new `rhoS`.
    8.  The change `DFr = MAXVAL(ABS(Fr-oldFr))` is checked against `scfc`.
    9.  If `DFr < scfc`, `Fr` is considered converged. The `pot_tol` is tightened.
    10. If not converged, `new_iteration` is set to `.TRUE.` to signal the calling SCF procedure to perform more cycles with the updated `Fr`.
  - **Arguments:**
    - `rho :: REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN)`: Current electron density.
    - `new_iteration :: LOGICAL, INTENT(INOUT)`: Output flag for SCF control.
    - `kinetic :: INTEGER, INTENT(INOUT)`: KEDF type selector.
    - `pot_tol :: REAL(KIND=DP), INTENT(INOUT)`: Potential tolerance.

- **`SUBROUTINE DecomposeDensityKEDF(rho,Fr,energy,potential,vaccutoff)`**:
  - **Description:** Calculates the WGCD kinetic energy and potential using the (converged) `Fr` matrix.
    It computes: `rho1 = rho * Fr`.
    `E_kin = (lambda*TF[rho] + mu*vW[rho]) - (lambda*TF[rho1] + mu*vW[rho1]) + (TF[rho1] + vW[rho1] + WGC[rho1; rhoS_from_rho1])`.
    The potential is derived accordingly. `WGC[rho1]` is calculated by calling `WGCPotentialPlus`.
  - **Arguments:**
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Total electron density.
    - `Fr :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: The converged scaling function `F(r)`.
    - `energy :: REAL(KIND=DP), INTENT(OUT)`: Calculated WGCD kinetic energy.
    - `potential :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT)`: Calculated WGCD kinetic potential.
    - `vaccutoff :: LOGICAL, INTENT(IN)`: Passed to `WGCPotentialPlus`.

# Important Variables/Constants

- **Module-Level Control & State:**
    - `WGCDiter :: INTEGER`: Iteration counter for the self-consistency of `Fr`.
    - `WGCDflag :: LOGICAL`: If `.TRUE.`, the WGCD iterative process for `Fr` is active.
    - `Fr(:,:,:) :: REAL(KIND=DP), ALLOCATABLE`: The spatially varying scaling function `F(r)`.
    - `oldFr(:,:,:) :: REAL(KIND=DP), ALLOCATABLE`: `Fr` from the previous `ComputeFrMatrix` iteration, used for convergence check.
- **Parameters for Scaling Function `F(rho/rho0 - shiftm)`:**
    - `t(151), scalefunt(151), scalefunDD(151) :: REAL(KIND=DP)`: Pre-tabulated x-values (`t`), function values (`scalefunt`), and second derivatives (`scalefunDD`) for spline interpolation of the function `f` that defines `Fr`.
    - `rhoc :: REAL(KIND=DP)`: A cutoff density used when `bvac` (vacuum treatment) is true, for determining `rho0` (default: 6.84e-3).
    - `shiftm :: REAL(KIND=DP)`: A shift parameter for the argument of the scaling function (default: 0.0).
- **Convergence Criterion for `Fr`:**
    - `scfc :: REAL(KIND=DP)`: Self-consistent field criterion for `Fr`. The `ComputeFrMatrix` loop terminates if `MAXVAL(ABS(Fr-oldFr)) < scfc` (default: 1e-4).
- **Imported KEDF Parameters:**
    - `lambda` (from `KEDF_TF`), `mu` (from `KEDF_VW`): Coefficients for the `T_1` functional.
    - `alpha`, `beta` (from `KEDF_WTkernel`): Density exponents, temporarily modified during `InitializeWGCD` and then restored.
- **Imported System Variables:**
    - `rho0`, `rhoS` (from `SYS`): Reference densities. `rho0` is recalculated iteratively in `ComputeFrMatrix`. `rhoS` is then set based on this new `rho0` for the WGC kernel generation.
    - `bvac` (from `SYS`): Flag for vacuum treatment.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP`, `PI`, `IMAG`.
- **`OutputFiles`, `Output`**: For `outputUnit`, `WrtOut`.
- **`MPI_Functions`**: For parallel communication (`rankGlobal`, `ReduceRealLevel1`, `MPI_ALLREDUCE`, `mpierr`) and timing (`StartClock`, `StopClock`).
- **`MathFunctions`**: For `stepfun`.
- **`MathSplines`**: For `spline_cubic_set` (though not explicitly called, `scalefunDD` implies it's used to set up the spline) and `spline_cubic_val` (for interpolating `Fr`).
- **`KEDF_TF`**: For the `lambda` parameter. The `DecomposeDensityKEDF` also implicitly uses TF energy and potential formulas.
- **`KEDF_VW`**: For the `mu` parameter. `DecomposeDensityKEDF` also implicitly uses vW energy and potential formulas.
- **`KEDF_WTkernel`**: For `keKernel` array, `fillWT` subroutine, and `alpha`, `beta` parameters. `InitializeWGCD` temporarily sets KEDF type to WT and uses its `alpha`/`beta`.
- **`KEDF_WGCkernel`**: For `FillWGC` subroutine (to regenerate WGC kernel with updated `rhoS`).
- **`KEDF_WGC`**: For `rhov` (vacuum density parameter for WGC's own cutoff, distinct from WGCD's `rhoc`) and `WGCPotentialPlus` (used in `DecomposeDensityKEDF` to calculate `T_2[rho1]`).
- **`SYS`**: For global reference densities `rho0`, `rhoS`, and `bvac` flag.
- **`Fourier`**: The `FFT` interface is used in `DecomposeDensityKEDF`.
- **`PlaneWave`**: `qTable` is used in `DecomposeDensityKEDF`.

**Workflow:**
1.  `InitializeWGCD` is called once at the start of an ion-optimization or single-point calculation if WGCD is selected. It sets the KEDF type temporarily to WT.
2.  The main SCF loop (e.g., in `RhoOptimizers`) proceeds. It will call `ComputeFrMatrix` in each outer SCF iteration (or when `new_iteration` is true).
3.  `ComputeFrMatrix` updates `Fr` and `rho0`. If `Fr` is not converged, it sets `new_iteration = .TRUE.`.
4.  The KEDF potential used by the SCF solver is calculated by `DecomposeDensityKEDF` using the current `Fr`.
5.  This cycle repeats until `Fr` converges. After `Fr` converges, `new_iteration` remains false, and the SCF proceeds normally with the fixed, converged `Fr` to optimize the electron density.
