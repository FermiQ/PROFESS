# Overview

The `KEDF_Q` module implements the LQ and HQ Kinetic Energy Density Functionals (KEDFs) developed by Jeng-Da Chai and John D. Weeks, as detailed in "Orbital-Free Density Functional Theory: Kinetic Potentials and Ab Initio Local Pseudopotentials," Phys. Rev. B 75, 205122 (2007).

These functionals are nonlocal and are designed to improve upon simpler local and semi-local KEDF approximations. They incorporate the Lindhard response function, `LindG(q/tkf, lambda, mu)`, where `tkf` is related to the Fermi wavevector of a reference density `rho0`, and `lambda` and `mu` are parameters often associated with the Thomas-Fermi and von Weizsacker functionals respectively. The LQ and HQ functionals differ in how they incorporate the density and its derivatives with the Lindhard function term.

# Key Components

- **`MODULE KEDF_Q`**: The main container for the LQ and HQ KEDF calculations.

- **`SUBROUTINE CalLHQ(potential, rho, calcEnergy, optSqrt, eTF, eVW, eQ, LQflag)`**:
  - **Description:** This is the primary driver routine for using either the LQ or HQ KEDF. It calculates the total KEDF potential and energy by summing three components:
    1.  Thomas-Fermi (TF) term (from `KEDF_TF` module).
    2.  von Weizsacker (vW) term (from `KEDF_VW` module).
    3.  The nonlocal Q-term, which is either `LQPotential`/`LQEnergy` or `HQPotential`/`HQEnergy` based on the `LQflag`.
  - **Arguments:**
    - `potential :: REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(OUT)`: Output KEDF potential.
    - `rho :: REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(IN)`: Input electron density.
    - `calcEnergy :: LOGICAL, INTENT(IN)`: If true, calculates energy components.
    - `optSqrt :: LOGICAL, INTENT(IN)`: If true, adjusts the potential for optimization w.r.t. `sqrt(rho)`.
    - `eTF, eVW, eQ :: REAL(KIND=DP), INTENT(OUT)`: Output energy components (TF, vW, Q-term).
    - `LQflag :: LOGICAL, INTENT(IN)`: If `.TRUE.`, the LQ functional is used; if `.FALSE.`, the HQ functional is used.

- **`FUNCTION LQEnergy(rho) RESULT(energy_val)`**:
  - **Description:** Calculates the specific nonlocal energy component of the LQ functional. It involves terms like `FFT(SQRT(rho))`, the Lindhard function `LindG`, and a real-space term `(3*rho - r . grad(rho))`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION LQPotential(rho) RESULT(potential_val)`**:
  - **Description:** Calculates the potential corresponding to the LQ nonlocal energy term. It involves `FFT(SQRT(rho))` convolved with the Lindhard function, then divided by `SQRT(rho)` in real space.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION HQEnergy(rho) RESULT(energy_val)`**:
  - **Description:** Calculates the specific nonlocal energy component of the HQ functional. It involves terms like `FFT(cTF*(5/3)*rho**(2/3))`, the Lindhard function, and a real-space term `(3*rho - r . grad(rho))`.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

- **`FUNCTION HQPotential(rho) RESULT(potential_val)`**:
  - **Description:** Calculates the potential corresponding to the HQ nonlocal energy term. It involves `FFT(cTF*(5/3)*rho**(2/3))` convolved with the Lindhard function.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`.

# Important Variables/Constants

- **Input/Control:**
    - `rho :: REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)`: The real-space electron density.
    - `LQflag :: LOGICAL, INTENT(IN)`: Selects between LQ and HQ functionals in `CalLHQ`.
- **Derived from `rho0` (via `SYS` module):**
    - `tkf :: REAL(KIND=DP)`: Calculated as `2.0 * (3.0 * PI**2 * rho0)**(1.0/3.0)`. This represents `2 * k_F(rho0)`, where `k_F` is the Fermi wavevector of the reference uniform electron gas density `rho0`.
- **Parameters from other KEDF modules:**
    - `lambda :: REAL(KIND=DP)` (from `KEDF_TF`): Coefficient associated with the Thomas-Fermi energy, used as a parameter in the `LindG` function.
    - `mu :: REAL(KIND=DP)` (from `KEDF_VW`): Coefficient associated with the von Weizsacker energy, used as a parameter in the `LindG` function.
    - `cTF :: REAL(KIND=DP)` (from `KEDF_TF`): The Thomas-Fermi constant `(3/10)*(3*PI^2)**(2/3)`.
- **Reciprocal Space Tools (from `PlaneWave` module):**
    - `qVectors :: REAL(KIND=DP), DIMENSION(k1G,k2G,k3G,3)`: Cartesian components of G-vectors.
    - `qTable :: REAL(KIND=DP), DIMENSION(k1G,k2G,k3G)`: Magnitudes of G-vectors, `|G|`.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`PlaneWave`**: For `qVectors` and `qTable` (G-vector information).
- **`Fourier`**: For the `FFT` interface, used to transform densities and potentials between real and reciprocal space.
- **`CONSTANTS`**: For `DP` (double precision), `PI`, and `IMAG` (imaginary unit).
- **`KEDF_TF`**: For `lambda` and `cTF` parameters, and for the base `TFPotential` and `TFEnergy` routines called by `CalLHQ`.
- **`KEDF_VW`**: For the `mu` parameter, and for the base `VWPotentialSqrtPlus` routine called by `CalLHQ`.
- **`SYS`**: For `rho0` (the reference density used to define `tkf`).
- **`CellInfo`**: For cell parameters (`cell%cellReal`), grid offsets (`n3Goff`), and global grid dimensions (`m3G`, `n1G`, `n2G`, `n3G`).
- **`MATHFUNCTIONS`**: For `LindG(q_reduced, lambda_param, mu_param)` which calculates the Lindhard response function or a similar kernel form.

**Workflow:**
When `CalLHQ` is called (e.g., from `CalPotPlus` if `kinetic` flag selects LQ or HQ KEDF):
1. It first calls `TFPotential` and `TFEnergy` for the Thomas-Fermi contribution.
2. Depending on `LQflag`, it calls either:
   - `LQPotential` and `LQEnergy` for the LQ nonlocal term.
   - `HQPotential` and `HQEnergy` for the HQ nonlocal term.
3. These nonlocal functions internally use `tkf` (derived from `rho0`), `lambda`, and `mu` as parameters for the `LindG` function, and perform FFTs for convolutions.
4. `CalLHQ` then calls `VWPotentialSqrtPlus` for the von Weizsacker contribution.
5. The potentials (and energies) from these three components are summed.
6. If `optSqrt` is true, the final potential is adjusted using the chain rule for `dE/d(sqrt(rho))`.
