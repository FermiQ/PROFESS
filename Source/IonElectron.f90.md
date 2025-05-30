# Overview

The `IonElectron` module calculates the contributions of the ion-electron interaction to the total energy, forces on ions, and the stress tensor. This interaction, often referred to as the external potential in Density Functional Theory (DFT), arises from the Coulomb attraction between the negatively charged electrons and the positively charged atomic nuclei (or ion cores, as represented by pseudopotentials).

The calculations are primarily performed in reciprocal space, leveraging the efficiency of Fast Fourier Transforms (FFTs) for convolutions that appear in real space. The module uses pseudopotentials to represent the ion cores, and structure factors to account for the positions of the ions in the unit cell.

# Key Components

- **`MODULE IonElectron`**: The main container for ion-electron interaction calculations.

- **`FUNCTION IonElectronEnergy(IonElectPotRecip, cellVol, rhoRecip) RESULT(energy_val)`**:
  - **Description:** Computes the ion-electron interaction energy using reciprocal space quantities. The formula is typically `E_ie = Omega * SUM_G { V_ei(G)^* * rho(G) }`, where `V_ei(G)` is the ion-electron potential in reciprocal space, `rho(G)` is the electron density in reciprocal space, and `Omega` is the cell volume. Special care is taken for the `G=0` term.
  - **Arguments:**
    - `IonElectPotRecip :: COMPLEX(kind=DP), DIMENSION(:,:,:), INTENT(IN)`: The ion-electron potential in reciprocal space.
    - `cellVol :: REAL(kind=DP), INTENT(IN)`: The volume of the simulation cell.
    - `rhoRecip :: COMPLEX(kind=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density in reciprocal space.
  - **Return Value:**
    - `energy_val :: REAL(kind=DP)`: The calculated ion-electron energy.

- **`FUNCTION IonElectronEnergyReal(rhoR_SI, pot) RESULT(energy_val)`**:
  - **Description:** Calculates the ion-electron energy in real space. The formula is `E_ie = SUM_r { rho(r) * V_ei(r) }`, where `rho(r)` is the real-space electron density and `V_ei(r)` is the real-space ion-electron potential.
  - **Arguments:**
    - `rhoR_SI :: REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density in real space.
    - `pot :: REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)`: The real-space ion-electron potential.
  - **Return Value:**
    - `energy_val :: REAL(kind=DP)`: The calculated ion-electron energy.

- **`FUNCTION IonElectronPotentialRecip(cellReal, ionTable, elementTable) RESULT(potential_recip)`**:
  - **Description:** Calculates the total ion-electron potential `V_ei(G)` in reciprocal space. This is constructed by summing the individual atomic pseudopotentials `V_ps(G; type)` multiplied by their respective structure factors `S_type(G) = SUM_I { exp(-iG.R_I_type) }`, and divided by the cell volume: `V_ei(G) = (1/Omega) * SUM_type { V_ps(G; type) * S_type(G) }`.
  - **Arguments:**
    - `cellReal :: REAL(KIND=DP), DIMENSION(3,3), INTENT(IN)`: Real-space lattice vectors.
    - `ionTable :: TYPE(ion), DIMENSION(:), INTENT(IN)`: Array of ion data (positions, types).
    - `elementTable :: TYPE(element), DIMENSION(:), INTENT(IN)`: Array of element data (pseudopotentials).
  - **Return Value:**
    - `potential_recip :: COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G)`: The ion-electron potential in reciprocal space.

- **`FUNCTION IonElectronForces(rhoRecip, ionTable, elementTable, cellReal) RESULT(forces_on_ions)`**:
  - **Description:** Computes the force on each ion due to the ion-electron interaction. The force on ion `I` is derived from the energy expression and is typically `F_I = SUM_G { iG * V_ps(G; type_I) * rho(G)^* * exp(iG.R_I) }`. (Note: The exact sign and complex conjugation details depend on FFT conventions and definitions of potential/structure factor).
  - **Arguments:**
    - `rhoRecip :: COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density in reciprocal space.
    - `ionTable :: TYPE(ion), DIMENSION(:), INTENT(IN)`: Array of ion data.
    - `elementTable :: TYPE(element), DIMENSION(:), INTENT(IN)`: Array of element data.
    - `cellReal :: REAL(KIND=DP), DIMENSION(3,3), INTENT(IN)`: Real-space lattice vectors.
  - **Return Value:**
    - `forces_on_ions :: REAL(KIND=DP), DIMENSION(SIZE(ionTable),3)`: Array containing the Cartesian components of the force on each ion.

- **`FUNCTION IonElectronStress(rhoRecip, ionTable, elementTable, cellReal, energy) RESULT(stress_tensor)`**:
  - **Description:** Calculates the contribution of the ion-electron interaction to the total stress tensor of the system. The formula involves a sum over reciprocal lattice vectors `G`, incorporating `rho(G)`, pseudopotential derivatives, structure factors, and outer products of `G`.
  - **Arguments:**
    - `rhoRecip :: COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: Spin-independent electron density in reciprocal space.
    - `ionTable :: TYPE(ion), DIMENSION(:), INTENT(IN)`: Array of ion data.
    - `elementTable :: TYPE(element), DIMENSION(:), INTENT(IN)`: Array of element data.
    - `cellReal :: REAL(KIND=DP), DIMENSION(3,3), INTENT(IN)`: Real-space lattice vectors.
    - `energy :: REAL(KIND=DP), INTENT(IN)`: The pre-calculated ion-electron energy (used for the trace part of the stress tensor).
  - **Return Value:**
    - `stress_tensor :: REAL(KIND=DP), DIMENSION(3,3)`: The 3x3 ion-electron stress tensor.

# Important Variables/Constants

- **`qTable :: REAL(KIND=DP), DIMENSION(k1G,k2G,k3G)`**: (Imported from `PlaneWave`) An array containing the magnitudes (`G`) of the reciprocal lattice vectors. Used for looking up pseudopotential values `V_ps(G)`.
- **`qMask :: LOGICAL, DIMENSION(k1G,k2G,k3G)`**: (Imported from `PlaneWave`) A mask indicating which G-vectors are active (e.g., within an energy cutoff).
- **`qVectors :: REAL(KIND=DP), DIMENSION(k1G,k2G,k3G,3)`**: (Imported from `PlaneWave`) An array containing the Cartesian components `(Gx, Gy, Gz)` for each G-vector.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision), `IMAG` (imaginary unit), `PI`. `auToGPa` is imported but not used in this specific file's functions.
- **`MathFunctions`**: For `Norm` (vector norm), `Vecmul` (matrix-vector product), `Inverse` (matrix inverse), `Volume` (cell volume).
- **`CellInfo`**: For data types (`ion`, `element`), the `cell` module variable (containing cell volume, etc.), and reciprocal space grid dimensions (`k1G, k2G, k3G, k3Goff, n3G, m3G`).
- **`LocalPseudoPot`**: For `PseudoPotLookup` (to get `V_ps(G)`) and `PseudoPotDiffLookup` (to get `dV_ps(G)/dG`, used in stress calculation).
- **`PlaneWave`**: For `qTable`, `qMask`, `qVectors`, and `CCStructureFactor` (calculates the complex conjugate of the ionic structure factor).
- **`OutputFiles`**: For `outputUnit` (Fortran I/O unit for standard output).
- **`IonElectronSpline`**: The module is used (`USE IonElectronSpline, ONLY: IonElectronPotRecipSpline`), but `IonElectronPotRecipSpline` is not directly called within the functions shown. This suggests that spline-based calculations might be handled by routines in `IonElectronSpline` itself, or this was a planned/alternative path.
- **`Fourier`**: (Implicit) The module deals with reciprocal space quantities (`rhoRecip`, `IonElectPotRecip`). These are typically obtained by FFTing their real-space counterparts using the `Fourier` module.

This module is central to computing the interaction between ions and electrons. Its outputs (energy, potential, forces, stress) are essential components for total energy evaluation, electronic structure optimization (via the potential), ionic force calculations, and cell optimization (via stress).
