# Overview

The `Ewald` module implements the Ewald summation technique to calculate the long-range ion-ion electrostatic interactions in a periodic system. This method is crucial for accurately determining energies, forces, and stresses in simulations of crystalline materials or other periodic structures.

The core idea of the Ewald sum is to split the slowly converging `1/r` Coulomb sum into two faster-converging parts:
1.  A **real-space sum**: This part handles short-range interactions directly. It involves a screened Coulomb potential, where a Gaussian charge distribution is added to each ion to screen its interaction, and an equal and opposite Gaussian is subtracted. The interaction of these screened charges converges rapidly.
2.  A **reciprocal-space sum**: This part accounts for the long-range interactions. It calculates the interaction of the smooth, subtracted Gaussian charge distributions via a sum over reciprocal lattice vectors. This also converges rapidly.

The module includes implementations for both the standard Ewald summation (scaling roughly as O(N^1.5) to O(N^2)) and the Particle Mesh Ewald (PME) method (scaling as O(N log N)), where "Spline" in function names usually indicates PME. PME uses FFTs and B-splines to efficiently calculate the reciprocal space sum.

The module handles the calculation of:
- Ewald energy
- Ewald forces on each ion
- Ewald contribution to the stress tensor

It also includes a setup routine (`EwaldSetup`) to determine optimal parameters for the summation (like the splitting parameter `eta` and cutoff radii) to achieve a desired accuracy.

**Key Assumptions:**
- The overall charge of the system is zero (charge neutrality).
- The system is fully periodic in three dimensions (no surfaces, or surface corrections are handled elsewhere).
- Ions are treated as point charges; this assumption can break down at very short interionic distances.

# Key Components

- **`MODULE Ewald`**: The main container for all Ewald-related subroutines and functions.

- **`SUBROUTINE EwaldSetup(cellReal, ionTable, elementTable, useSpline, sizeX, sizeY, sizeZ)`**:
  - **Description:** Calculates and sets up the optimal Ewald parameters: `eta` (the Ewald splitting parameter), `realCutoff` (real-space cutoff), and `recipCutoff` (reciprocal-space cutoff). These parameters are determined based on the desired `errorTolerance`, cell dimensions, and whether PME is used. For PME, `sizeX, sizeY, sizeZ` (dimensions of the FFT grid) are also relevant.
  - **Key actions:** Initializes module-level variables like `cellVol`, `sumCharge`, `sumChargeSquared`, `cellRecip`, `eta`, `realCutoff`, `recipCutoff`, `nArecip`, `nBrecip`, `nCrecip`, and their PME counterparts (`nArecipSpl`, etc.).

- **`FUNCTION EwaldEnergy(cellReal, ionTable, elementTable, useSpline) RESULT(totalEwaldEnergy)`**:
  - **Description:** Computes the total ion-ion Ewald interaction energy. It sums contributions from real space, reciprocal space, a self-interaction correction, and an average potential term.
  - **Contains Internal Functions:**
    - `FUNCTION EwaldRealSpace(eta, realCutoff)`: Calculates the real-space sum part of the energy.
    - `FUNCTION EwaldRealSpaceSelfTerm(eta)`: Calculates the self-energy correction (interaction of an ion's Gaussian charge with itself).
    - `FUNCTION EwaldRecipSpace(eta, recipCutoff)`: Calculates the reciprocal-space sum part of the energy using direct summation over G-vectors.
    - `FUNCTION EwaldRecipSpaceSpline(eta, recipCutoff)`: Calculates the reciprocal-space sum using the PME method.
    - `FUNCTION EwaldAveragePotentialTerm(eta)`: Calculates a correction term related to the average potential of the unit cell or the neutralizing background.
  - **Output:** `ionIonEnergy` (module variable) is updated with the total Ewald energy.

- **`FUNCTION EwaldForces(cellReal, ionTable, elementTable, useSpline) RESULT(forcesMatrix)`**:
  - **Description:** Computes the forces on each ion due to ion-ion Ewald interactions.
  - **Contains Internal Functions/Subroutines:**
    - `FUNCTION EwaldRealSpaceForces(iCoord, iCharge, eta)`: Calculates the real-space force contribution on a specific ion `i`.
    - `SUBROUTINE EwaldRecipSpaceForces(eta, recipCutoff, EwaldForces)`: Calculates the reciprocal-space force contributions for all ions. (Note: This was changed from a FUNCTION to a SUBROUTINE that updates `EwaldForces`.)
    - `SUBROUTINE EwaldRecipSpaceForcesSpline(eta, recipCutoff, EwaldForces)`: Calculates reciprocal-space forces using PME.

- **`FUNCTION EwaldStress(cellReal, ionTable, elementTable, useSpline) RESULT(stressTensor)`**:
  - **Description:** Computes the ion-ion Ewald contribution to the system's stress tensor.
  - **Contains Internal Functions:**
    - `FUNCTION EwaldRealSpaceStress(eta)`: Calculates the real-space contribution to the stress tensor.
    - `FUNCTION EwaldRecipSpaceStress(eta, recipCutoff)`: Calculates the reciprocal-space contribution to the stress tensor.
    - `FUNCTION EwaldRecipSpaceStressSpline(eta, recipCutoff)`: Calculates the reciprocal-space stress using PME.
    - `FUNCTION EwaldAveragePotentialTermStress(eta)`: Calculates the stress contribution from the average potential term.

# Important Variables/Constants

**User-Configurable Parameters (Module-Level):**
- `errorTolerance :: REAL(KIND=DP)`: Desired precision for Ewald sums (default: `1.e-10 / hartreeToeV`, i.e., ~1e-10 Ha).
- `maxRealCutoff :: REAL(KIND=DP)`: Maximum allowed real-space cutoff distance (default: `12.0 / BOHR`, i.e., 12 Bohr).
- `etaIncrement :: REAL(KIND=DP)`: Step size for optimizing `eta` (default: `0.01 * BOHR`).
- `recipCutoffIncrement :: REAL(KIND=DP)`: Step size for optimizing `recipCutoff` (default: `0.01 * BOHR`).
- `realCutoffIncrement :: REAL(KIND=DP)`: Step size for optimizing `realCutoff` (default: `0.01 * BOHR`).

**Calculated Ewald Parameters & Shared Data (Module-Level, Private):**
- `eta :: REAL(KIND=DP)`: The Ewald splitting parameter, balancing real and reciprocal space work.
- `realCutoff :: REAL(KIND=DP)`: Optimized real-space cutoff radius.
- `recipCutoff :: REAL(KIND=DP)`: Optimized reciprocal-space cutoff radius.
- `cellVol :: REAL(KIND=DP)`: Volume of the simulation cell.
- `sumCharge :: REAL(KIND=DP)`: Sum of all ionic charges in the cell.
- `sumChargeSquared :: REAL(KIND=DP)`: Sum of squared ionic charges.
- `sqrtPi :: REAL(KIND=DP)`: Square root of pi.
- `cellRecip :: REAL(KIND=DP), DIMENSION(3,3)`: Reciprocal lattice vectors.
- `nArecip, nBrecip, nCrecip :: INTEGER`: Number of reciprocal lattice vector repetitions for standard Ewald sum.
- `nArecipSpl, nBrecipSpl, nCrecipSpl :: INTEGER`: Number of grid points for PME FFT grid (taken from `sizeX, sizeY, sizeZ` input to `EwaldSetup`).
- `numIon :: INTEGER`: Total number of ions.
- `numIonType :: INTEGER`: Number of different ion types.
- `ionScreen :: LOGICAL`: Flag for using efficient ion screening (default `.TRUE.`).
- `ionIonEnergy :: REAL(KIND=DP)`: Stores the calculated total ion-ion Ewald energy.

**Temporary Arrays (Module-Level, Allocatable, Private - for PME):**
- `tempRA`: General-purpose real-space work array.
- `qIonTable`: Array `Q` from Essmann et al. [5], charge assigned to grid points.
- `CCStructure`: Complex conjugate of the structure factor `S(k)`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

This module relies on several other modules:

- **`Constants`**: For `DP`, `PI`, `hartreeToeV`, `BOHR`.
- **`CellInfo`**: For data types (`ion`, `element`), FFT grid dimensions (`k1G`, `k2G`, `k3G`, `k3Goff`, `n1G`, `n2G`, `n3G`, `n3Goff`, `m1G`, `m2G`, `m3G`, `m123G`).
- **`Fourier_NEW`**: For FFT routines (`FFT_NEW`, `FFT_STD_STATE`).
- **`MathFunctions`**: For vector/matrix operations (`Vecmul`, `Norm`, `Inverse`, `Volume`) and `Metric3d`.
- **`Timer`**: For performance timing (`TimerStart`, `TimerStop`, `stopwatch` type).
- **`OUTPUT`**: For writing output (`WrtOut`, `PrintEwald`) and logging (`Title`, `StartClock`, `StopClock`).
- **`OutputFiles`**: For `outputUnit`.
- **`MPI_Functions`**: For parallel communication (`MPI_ALLREDUCE`, `MPI_COMM_RANK`, `MPI_COMM_SIZE`, etc.) if `__USE_PARALLEL` is defined.
- **`SYS`**: For parallel-specific variables like `numIonLoc`, `numIonInit` if `__USE_PARALLEL` is defined.
- **`CBSpline`**: For PME calculations, specifically `BSplineProduct` and `splineOrder`.
- **`IonElectronSpline`**: For PME, `FillQIonTable` (to map ionic charges to grid) and `CalculateSplineForces`.
- **`PlaneWave`**: For `qVectors` (reciprocal space vectors corresponding to grid points).

The typical workflow involving this module is:
1.  Call `EwaldSetup` once the cell parameters and ion types are known, or if the cell changes significantly. This optimizes `eta`, `realCutoff`, and `recipCutoff`.
2.  Call `EwaldEnergy` to get the ion-ion energy.
3.  Call `EwaldForces` to get the ion-ion forces.
4.  Call `EwaldStress` to get the ion-ion contribution to the stress tensor.

The `useSpline` flag determines whether standard Ewald or PME methods are used for the reciprocal space parts. Parallel execution is handled via MPI calls bracketed by `#ifdef __USE_PARALLEL`.
