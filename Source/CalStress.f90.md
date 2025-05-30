# Overview

The `CalStress` module is dedicated to calculating the symmetric 3x3 stress tensor for the simulated system. It aggregates contributions from various physical interactions, including:
- Kinetic Energy Density Functional (KEDF) stress
- Hartree stress (from electron-electron Coulomb interaction)
- Exchange-correlation (XC) stress (supporting LDA or PBE functionals)
- Ion-electron stress (with an option for spline-based approximation)
- Ion-ion stress (typically from Ewald summation, also with a spline option)

The module acts as a driver, invoking specific stress calculation routines from other specialized modules.

# Key Components

- **`MODULE CalStress`**:
  - The main container for the stress calculation subroutine.

- **`SUBROUTINE CalculateStress(rhoR, energy, stress)`**:
  - **Description:** This is the core routine that computes the total stress tensor. It assumes that the electronic energy components (passed via the `energy` array) have already been calculated for the given electron density `rhoR`. The routine sums up stress contributions from KEDF, ion-electron interaction, Hartree potential, exchange-correlation, and ion-ion interaction. Finally, it ensures the resulting stress tensor is symmetric.
  - **Arguments:**
    - `rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN), TARGET`: The input electron density in real space, potentially spin-dependent.
    - `energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`: An array containing pre-calculated energy components of the system. Specific indices are used for different KEDF terms (e.g., `energy(7)` for TF, `energy(9)` for WT/WGC nonlocal parts), Hartree energy (`energy(4)`), XC energy (`energy(5)`), and ion-electron energy (`energy(3)`).
    - `stress :: REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT)`: The output 3x3 stress tensor.

# Important Variables/Constants

- **`rhoR :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN), TARGET`**: The real-space electron density, indexed by grid points and spin.
- **`energy :: REAL(KIND=DP), DIMENSION(:), INTENT(IN)`**: Array holding various energy components. The stress calculation for some terms (e.g., Thomas-Fermi) directly uses these pre-calculated energies.
    - `energy(3)`: Ion-electron energy.
    - `energy(4)`: Hartree energy.
    - `energy(5)`: Exchange-correlation energy.
    - `energy(7)`: Thomas-Fermi KEDF component of energy.
    - `energy(9)`: Wang-Teter or Wang-Govind-Carter non-local KEDF energy component.
- **`stress :: REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT)`**: The resulting 3x3 stress tensor.
- **`rhoRecip_SI :: COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G,SIZE(rhoR,4))`**: Reciprocal space representation of the electron density.
- **`kinetic :: INTEGER`**: (from `CellInfo` module) Integer flag selecting the KEDF to be used for the kinetic stress contribution.
- **`exchangeCorrelation :: INTEGER`**: (from `XC_LDA` module) Integer flag selecting the XC functional (LDA, PBE) for the XC stress contribution.
- **`iiSpline :: LOGICAL`**: (from `IonElectronSpline` module) Boolean flag controlling whether a spline approximation is used for ion-ion interactions within the Ewald stress calculation.
- **`ieSpline :: LOGICAL`**: (from `IonElectronSpline` module) Boolean flag controlling whether a spline approximation is used for the ion-electron stress calculation.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

This module collaborates with numerous other modules to calculate the total stress:

**Core & Utility Modules:**
- `Constants`: For `DP` (double precision) and `auToGPa` (conversion factor).
- `CellInfo`: For cell parameters (`cell`), number of spins (`numSpin`), KEDF selector (`kinetic`), and reciprocal space grid info (`k1G`, `k2G`, `k3G`, `m123G`).
- `Output`: For utility subroutines like `QUIT`, `TITLE`, `StartClock`, `StopClock`. (Note: `TITLE`, `StartClock`, `StopClock` might be part of `MPI_Functions` or a similar utility if not directly in `Output`).
- `MPI_Functions`: For `ReduceRealLevel1` (MPI reduction for summing stress contributions in parallel).
- `OutputFiles`: General use, no specific entities mentioned in the direct context of `CalculateStress`.
- `Fourier_NEW`: For `FFT_STD_STATE`, `FFT_NEW` (FFT routines).
- `Fourier`: For `FFT` (another FFT routine).

**Physics Modules (Stress Components):**
- `IonElectronSpline`: For boolean flags `iiSpline`, `ieSpline` and the `IonElectronStressSpline` subroutine.
- `IonElectron`: For the `IonElectronStress` subroutine.
- `Hartree`: For `JStress` (Hartree stress).
- `Ewald`: For `EwaldStress` (ion-ion Ewald stress).

**Exchange-Correlation Stress Modules:**
- `XC_PBE`: For `PBEStress` (PBE GGA stress).
- `XC_LDA`: For `exchangeCorrelation` selector, `LDAStress`, `LSDAStress` (LDA and LSDA stress).

**KEDF Stress Modules:**
- `KEDF_TF`: For `TFStress` (Thomas-Fermi stress).
- `KEDF_VW`: For `VWStress` (von Weizsaecker stress).
- `KEDF_WT`: For `WTStress` (Wang-Teter stress).
- `KEDF_WGC`: For `WGCStress` (Wang-Govind-Carter stress).
- `KEDF_HC10`: For `IntStress` (Huang-Carter KEDF stress).

The `CalculateStress` subroutine sequentially calls these specialized stress functions, accumulates their contributions, performs an MPI reduction if in parallel, adds the Ewald stress, and finally symmetrizes the total stress tensor. The selection of KEDF and XC functional is driven by the integer flags `kinetic` and `exchangeCorrelation`.
