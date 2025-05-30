# Overview

The `KEDF_DenDec` module implements the density decomposition method for orbital-free DFT, as proposed by Huang and Carter in Phys. Rev. B 85, 045126 (2012). This approach is particularly aimed at improving the accuracy of KEDFs for systems like transition metals.

The core concept is to decompose the total electron density `rho_total(r)` into two parts:
1.  `rho_core(r)`: A localized, core-like density component that is typically constructed from atomic core densities and is assumed to follow the nuclei rigidly (though Pulay forces due to its ion dependence are considered).
2.  `rho_delocalized(r)`: The remaining delocalized or "valence" electron density, which is the primary variable optimized in the OFDFT calculation.

The kinetic energy is then expressed not just as `T[rho_total]` but often involves terms like `T[rho_delocalized]` plus correction terms that depend on both `rho_core` and `rho_delocalized` (or `rho_total`). This module provides routines to manage these densities and calculate the specific potential and force contributions arising from this decomposition.

# Key Components

- **`MODULE KEDF_DenDec`**: The main module for density decomposition functionalities.

- **`SUBROUTINE SetupCoreDensity`**:
  - **Description:** Initializes the total core density `core_den` for the system. It calls `MakeCoreDensity` to construct `core_den` from atomic core density files. Then, it adjusts the initial delocalized density (`rhoR` from `Sys` module) by subtracting the integrated charge of `core_den` and rescaling `rhoR` to ensure the total electron number matches the sum of ionic charges.
  - **Triggered by:** `do_den_dec == 1`.

- **`SUBROUTINE AllocateCoreDensity`**:
  - **Description:** Allocates module-level arrays: `core_den` (real space total core density), `potDD` (real space potential for Pulay forces), and `core_den_recip` (reciprocal space total core density). Sizes are determined by the simulation grid dimensions.
  - **Triggered by:** `do_den_dec == 1`.

- **`SUBROUTINE MakeCoreDensity`**:
  - **Description:** Constructs the total core density `core_den`. For each ion type with a specified `atomicCoreFile`:
    1. Reads 1D atomic core density data (typically `rho_atomic_core(q)`) from the file.
    2. Uses spline interpolation (`pchez`, `pchev` from `NMS` module) to evaluate `rho_atomic_core(G)` for each G-vector magnitude.
    3. Multiplies by the ionic structure factor `S_type(G)` (from `CCStructureFactor`) to get the contribution to `core_den_recip(G)`.
    4. After summing contributions from all ion types, `core_den_recip` is inverse FFTed to obtain the real-space `core_den`.
    5. Ensures `core_den` is non-negative.

- **`SUBROUTINE SubstrateCoreDensity(rhoR)`**: (Intended: SubtractCoreDensity)
  - **Description:** Subtracts the constructed `core_den` from the provided `rhoR` (total density) to obtain the delocalized density. `rhoR_new = rhoR_old - core_den`.
  - **Arguments:** `rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`.

- **`SUBROUTINE AddCoreDensity(rhoR)`**:
  - **Description:** Adds the `core_den` to the provided `rhoR` (delocalized density) to obtain the total density. `rhoR_new = rhoR_old + core_den`.
  - **Arguments:** `rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`.

- **`SUBROUTINE PulayForceDenDec(forces)`**:
  - **Description:** Calculates the Pulay force contribution that arises because `core_den` depends on ion positions. The force is essentially `integral{ potDD(r) * grad_R[core_den(r-R)] dr }`. `potDD` is the KEDF potential derivative w.r.t. `rho_delocalized`, and `grad_R[core_den]` is computed via FFTs of `iG * core_den_recip(G)`. (This seems to be an older or alternative version of `ForceDenDec`).
  - **Arguments:** `forces :: REAL(KIND=DP), DIMENSION(:,:), INTENT(INOUT)`.

- **`SUBROUTINE ForceDenDec(forces)`**:
  - **Description:** Calculates the Pulay-like forces due to the ion-dependent frozen core density. For each ion `I` of type `T`, it computes `F_I = SUM_G { CONJG(potDD(G)) * iG * rho_atomic_core_T(G) * exp(-iG.R_I) }`. This sums the force contribution from each G-vector component.
  - **Arguments:** `forces :: REAL(KIND=DP), DIMENSION(:,:), INTENT(INOUT)`.

- **`SUBROUTINE PotDenDec(rhoReal_SI, del_rho, potential, locETable, calcEnergy)`**:
  - **Description:** Calculates the specific KEDF potential term arising from the density decomposition. This usually involves a correction term added to the main KEDF potential. The correction is typically of the form `aTF * (V_TF[rho_total] - V_TF[rho_delocalized]) + bVW * (V_vW[rho_total] - V_vW[rho_delocalized])`, where `rho_total` is `rhoReal_SI` and `rho_delocalized` is `del_rho`. Adjustments are made if the optimization variable is `sqrt(rho_delocalized)`.
  - **Arguments:**
    - `rhoReal_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Total spin-independent density.
    - `del_rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Delocalized spin-independent density.
    - `potential :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`: KEDF potential array to be modified.
    - `locETable :: REAL(KIND=DP), DIMENSION(9), INTENT(INOUT)`: Array to store energy components.
    - `calcEnergy :: LOGICAL, INTENT(IN)`: Flag to compute energy.

# Important Variables/Constants

- **`do_den_dec :: INTEGER`**: Module-level flag. If 1, density decomposition is performed; if 0, this module's routines are largely bypassed.
- **`ncore :: INTEGER`**: Number of q-points in the tabulated 1D atomic core density files.
- **`recip_core(:,:,:) :: REAL(KIND=DP), ALLOCATABLE`**: Stores the 1D `q` and `rho_atomic_core(q)` data for each ion type read from `atomicCoreFile`. Shape: `(numIonTypes, ncore, 2)`.
- **`core_den(:,:,:) :: REAL(KIND=DP), ALLOCATABLE`**: The total constructed core electron density in real space on the simulation grid.
- **`core_den_recip(:,:,:) :: COMPLEX(KIND=DP), ALLOCATABLE`**: The total constructed core electron density in reciprocal space.
- **`potDD(:,:,:) :: REAL(KIND=DP), ALLOCATABLE`**: The real-space potential used to calculate Pulay forces. This is `delta E_correction / delta rho_delocalized`, where `E_correction` is the energy term like `aTF * (E_TF[rho_total] - E_TF[rho_delocalized]) + ...`.
- **`aTF :: REAL(KIND=DP)`**: Coefficient for the Thomas-Fermi part of the KEDF correction term (default 0.8).
- **`bVW :: REAL(KIND=DP)`**: Coefficient for the von Weizsaecker part of the KEDF correction term (default 0.2).
- **`atomicCoreFile(:) :: CHARACTER(LEN=100), ALLOCATABLE`**: Array storing the filenames for the 1D atomic core densities for each element type.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`, `PI`, `imag`, `hartreeToeV`, `bohr`.
- **`MATHFUNCTIONS`**: For `Volume`, `Inverse`, `Vecmul`, `Norm`.
- **`OutputFiles`**: For `outputUnit`.
- **`MPI_Functions`**: For MPI communication (`ReduceRealLevel1`, `ReduceReal_array`) and utility functions (`TITLE`, `StartClock`, `StopClock`).
- **`CellInfo`**: For cell parameters (`cell`), grid dimensions (`m123G`, `k1G`, `k2G`, `k3G`), and `numSpin`.
- **`Sys`**: For `rhoR` (the delocalized electron density).
- **`PlaneWave`**: For `qMask`, `qVectors` (G-vector information), and `CCStructureFactor` (ionic structure factor).
- **`FOURIER`**: For the `FFT` interface, used to transform densities and potentials between real and reciprocal space.
- **`NMS`**: For spline interpolation routines `pchez` (setup) and `pchev` (evaluation), used in `MakeCoreDensity` to evaluate atomic core densities at arbitrary q-points.
- **`KEDF_TF`**: For `TFPotential` and `TFEnergy` (used in `PotDenDec`).
- **`KEDF_VW`**: For `VWPotentialSqrtPlus` (used in `PotDenDec`).

When `do_den_dec` is active, this module significantly alters the KEDF calculation:
1.  `SetupCoreDensity` is called during initialization to construct `core_den` and adjust the initial `rhoR`.
2.  During self-consistency or optimization, the KEDF potential calculation (e.g., in `CalPotPlus`) will call `PotDenDec` to add the decomposition-specific correction potential.
3.  Force calculations (e.g., in `CalForces`) will call `ForceDenDec` (or `PulayForceDenDec`) to add the forces arising from the ion-dependent `core_den`.
4.  Routines like `AddCoreDensity` and `SubstrateCoreDensity` (SubtractCoreDensity) are used to convert between total density and delocalized density when interacting with parts of the code that expect one or the other.
