# Overview

The `KEDF_TF` module provides functions to calculate the Thomas-Fermi (TF) component of the kinetic energy, its corresponding potential, and its contribution to the stress tensor. The Thomas-Fermi model is the simplest local density approximation (LDA) for the kinetic energy of a non-interacting electron gas.

The TF kinetic energy is given by the formula:
`E_TF = lambda * cTF * integral( rho(r)^(5/3) dr )`
where:
- `rho(r)` is the electron density.
- `cTF` is the Thomas-Fermi constant, equal to `(3/10) * (3 * pi^2)^(2/3)`.
- `lambda` is a scaling parameter, often 1.0 in the pure TF model, but can be adjusted in more complex functionals (e.g., to represent only a fraction of the TF energy).

The corresponding potential is `V_TF(r) = dE_TF / drho(r) = lambda * (5/3) * cTF * rho(r)^(2/3)`.

# Key Components

- **`MODULE KEDF_TF`**: The main container for TF KEDF calculations.

- **`SUBROUTINE CalTF(potential, rho, calcEnergy, eTF)`**:
  - **Description:** Calculates the Thomas-Fermi potential and, optionally, the energy. It handles both spin-unpolarized (`numSpin=1`) and spin-polarized (`numSpin=2`) cases. For spin-polarized cases, the TF functional is typically applied to `2*rho_spin` for each spin channel, and the energies are summed (multiplied by 0.5).
  - **Arguments:**
    - `potential :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT)`: Output TF potential (added to the input `potential`).
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: Input electron density (4th dimension for spin).
    - `calcEnergy :: LOGICAL, INTENT(IN)`: If `.TRUE.`, calculates the TF energy.
    - `eTF :: REAL(KIND=DP), INTENT(OUT)`: Calculated TF energy.

- **`FUNCTION TFEnergy(rhoR_SI, weightingfunc, IfCalForce) RESULT(energy_val)`**:
  - **Description:** Calculates the Thomas-Fermi kinetic energy for a given spin-independent real-space density `rhoR_SI`. An optional `weightingfunc` can be provided.
  - **Arguments:**
    - `rhoR_SI :: REAL(KIND=DP),DIMENSION(:,:,:),INTENT(IN)`: Spin-independent electron density.
    - `weightingfunc :: REAL(KIND=DP),DIMENSION(:,:,:),INTENT(IN),OPTIONAL`: An optional weighting function.
    - `IfCalForce :: LOGICAL, INTENT(IN),OPTIONAL`: If present and true, `weightingfunc` is applied.
  - **Return Value:** `energy_val :: REAL(KIND=DP)`: The TF energy.

- **`FUNCTION TFPointEnergy(rhoR_SI, i, j, k) RESULT(point_energy_val)`**:
  - **Description:** Calculates the Thomas-Fermi kinetic energy density (`lambda * cTF * rho(i,j,k)^(5/3)`) at a single grid point `(i,j,k)`.
  - **Arguments:**
    - `rhoR_SI :: REAL(kind=DP), DIMENSION(0:,0:,0:), INTENT(IN)`: Spin-independent electron density.
    - `i, j, k :: INTEGER, INTENT(IN)`: Grid indices.
  - **Return Value:** `point_energy_val :: REAL(kind=DP)`.

- **`FUNCTION TFPotential(rhoR_SI) RESULT(potential_val)`**:
  - **Description:** Calculates the Thomas-Fermi potential for a spin-independent density `rhoR_SI`. `V_TF = (5/3 * cTF * lambda) * rhoR_SI**(2/3)`.
  - **Arguments:** `rhoR_SI :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.
  - **Return Value:** `potential_val :: REAL(KIND=DP), DIMENSION(SIZE(rhoR_SI,1), ...)`: The TF potential.

- **`FUNCTION TFPotentialSqrt(rhoR_SI, weightingfunc, IfCalForce) RESULT(potential_sqrt_val)`**:
  - **Description:** Calculates the Thomas-Fermi potential derivative with respect to `sqrt(rhoR_SI)`. This is `dE_TF / d(sqrt(rho)) = (dE_TF/drho) * (drho/d(sqrt(rho))) = V_TF * 2*sqrt(rho) = (10/3 * cTF * lambda) * rhoR_SI**(7/6)`.
  - **Arguments:** Same as `TFEnergy`.
  - **Return Value:** `potential_sqrt_val :: REAL(KIND=DP), DIMENSION(SIZE(rhoR_SI,1), ...)`: The TF potential w.r.t `sqrt(rho)`.

- **`FUNCTION TFStress(cellVol, energy) RESULT(stress_tensor)`**:
  - **Description:** Calculates the Thomas-Fermi contribution to the stress tensor. For TF, the stress is isotropic (diagonal) and given by `Stress_diag = - (2/3) * E_TF / cellVol`.
  - **Arguments:**
    - `cellVol :: REAL(KIND=DP), INTENT(IN)`: Volume of the cell.
    - `energy :: REAL(KIND=DP), INTENT(IN)`: The pre-calculated TF energy.
  - **Return Value:** `stress_tensor :: REAL(KIND=DP), DIMENSION(3,3)`: The TF stress tensor.

# Important Variables/Constants

- **Module-Level Parameters:**
    - **`lambda :: REAL(KIND=DP)`**: A scaling factor for the TF contribution. Default value is -100.0, which is unphysical and indicates that this parameter **must be initialized** to a sensible value (e.g., 1.0 for pure TF, or other values like 5/6 when TF is part of a GGA KEDF) elsewhere in the code, likely during the setup phase (e.g., in `SetupKEDF`).
    - **`cTF :: REAL(KIND=DP), PARAMETER`**: The Thomas-Fermi constant, `(3/10)*(3*pi^2)^(2/3)`, with a value of `2.87123400018819`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` (double precision kind parameter).
- **`PlaneWave`**: Imports `qTable`, `qVectors`, `qMask`. However, these reciprocal space quantities are not directly used by the TF functional itself, as TF is a purely local functional depending only on `rho(r)` at each point. Their import might be a remnant or for consistency with other KEDF modules.
- **`CellInfo`**: For `cell` (cell information), `numSpin` (used in `CalTF`), and grid parameters `n3G`, `m3G` (used in `TFStress` for parallel scaling).

The TF KEDF is often used as a base component in more advanced KEDFs (like GGAs or nonlocal functionals). Routines from this module (e.g., `TFPotential`, `TFEnergy`) are called by other KEDF modules (e.g., `KEDF_Q`, `KEDF_EvW`) to get the TF part of the energy/potential. The `lambda` parameter allows these other functionals to use a specific fraction of the TF energy/potential.
