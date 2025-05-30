# Overview

The `RhoDirBFGS` module implements the Limited-memory BFGS (L-BFGS) algorithm specifically for determining the search direction during the optimization of the electron density `phi` (where `phi` is typically `sqrt(rho)` in many orbital-free DFT approaches). L-BFGS is a quasi-Newton method that iteratively builds an approximation to the inverse Hessian matrix using information from a limited number (`MBFGS`) of previous steps. This avoids the need to store or invert a full Hessian, making it suitable for large-scale problems.

The module provides routines to:
1.  Initialize and manage the storage for the L-BFGS history vectors (`s_k = phi_{k+1} - phi_k` and `y_k = gradient_{k+1} - gradient_k`).
2.  Calculate a new search direction `dirNext = -H_k * gradient_k` using the L-BFGS two-loop recursion, where `H_k` is the approximate inverse Hessian.
3.  Update the L-BFGS history vectors after a successful step along the search direction.
4.  Calculate the chemical potential, which is needed to ensure the gradient `dL/dPhi` (where `L = E - mu * N_electrons`) is correctly defined for the constrained optimization of `phi` under electron number conservation.

**Reference:**
The L-BFGS algorithm implemented here is likely based on standard descriptions, such as the one found in Nocedal and Wright, "Numerical Optimization," or earlier papers like Liu and Nocedal, "On the limited memory BFGS method for large scale optimization," Math. Programming 45, 503-528 (1989). The comment mentions Math. of Comp. 35, 151, 1980, which might refer to earlier BFGS work.

# Key Components

- **`MODULE RhoDirBFGS`**: The main container module.

- **Module-Level Parameters & Variables:**
    - `MBFGS :: INTEGER`: The number of previous steps (pairs of `s` and `y` vectors) stored to approximate the Hessian (default: 5).
    - `BFGSPT :: INTEGER`: A pointer or index indicating the current slot in the circular buffers used to store the `s` and `y` vectors.
    - `BFGSiter :: INTEGER`: Tracks the current BFGS iteration number.
    - `BFGSrho(MBFGS, nspin) :: REAL(KIND=DP), ALLOCATABLE`: Stores the scalar values `rho_k = 1 / (y_k^T s_k)`.
    - `BFGSalpha(MBFGS, nspin) :: REAL(KIND=DP), ALLOCATABLE`: Stores the temporary scalar values `alpha_i` computed during the first loop of the L-BFGS two-loop recursion.
    - `BFGSs(MBFGS, dimX,dimY,dimZ,nspin) :: REAL(KIND=DP), ALLOCATABLE`: Stores the `MBFGS` most recent `s_k = phi_{k+1} - phi_k` vectors.
    - `BFGSy(MBFGS, dimX,dimY,dimZ,nspin) :: REAL(KIND=DP), ALLOCATABLE`: Stores the `MBFGS` most recent `y_k = (dL/dPhi)_{k+1} - (dL/dPhi)_k` vectors.

- **`SUBROUTINE InitializeBFGS(rho)`**:
  - **Description:** Allocates and initializes the module-level arrays (`BFGSrho`, `BFGSalpha`, `BFGSs`, `BFGSy`) based on the dimensions of the input density `rho` (to get grid size and number of spin channels `nspin`) and the parameter `MBFGS`. It also resets `BFGSPT` and `BFGSiter` to zero.
  - **Arguments:** `rho :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`.

- **`SUBROUTINE CleanBFGS()`**:
  - **Description:** Deallocates all allocatable module-level arrays used for L-BFGS history.

- **`SUBROUTINE BFGSDirection(energy, phi, dLdPhi, dLdPhi_old, dV, dirNext)`**:
  - **Description:** Calculates the new L-BFGS search direction `dirNext`.
    - If it's the first iteration (`BFGSiter == 1`), it returns the steepest descent direction: `dirNext = -dLdPhi`.
    - Otherwise, it implements the L-BFGS two-loop recursion:
        1.  First loop (backward): Computes `alpha_i = rho_i * s_i^T * q` and updates `q = q - alpha_i * y_i`.
        2.  Scales the result `q` by an initial Hessian approximation `H_0` (here, `H_0` is taken as a diagonal matrix `Diag = (s_k^T y_k) / (y_k^T y_k) * I`).
        3.  Second loop (forward): Computes `beta_i = rho_i * y_i^T * r` and updates `r = r + s_i * (alpha_i - beta_i)`.
    The final `r` (stored in `temp` and then `dirNext`) is the L-BFGS direction `-H_k * dLdPhi`.
  - **Arguments:**
    - `energy :: REAL(KIND=DP), INTENT(IN)`: Current total energy.
    - `phi :: REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:)`: Current `phi` variables.
    - `dLdPhi :: REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:)`: Current gradient `dL/dPhi`.
    - `dLdPhi_old :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT)`: Previous gradient (updated to current `dLdPhi` at the end).
    - `dV :: REAL(KIND=DP), INTENT(IN)`: Real-space grid volume element (used by `InnerProd`).
    - `dirNext :: REAL(KIND=DP),DIMENSION(SIZE(phi)...), INTENT(OUT)`: The calculated search direction.

- **`SUBROUTINE UpdateBFGS(phi, tempPhi, dEdPhi, dLdPhi, totEleNum, dV)`**:
  - **Description:** Updates the L-BFGS history vectors (`BFGSs`, `BFGSy`) and `BFGSrho` after a successful line search step.
    1.  Calls `ChemicalPotential` to calculate the current chemical potential `mu` based on `dEdPhi` (gradient of energy w.r.t. `phi`) and `tempPhi` (the new `phi` after line search).
    2.  Computes the updated gradient `tempdLdPhi = dEdPhi - 2*mu*tempPhi`.
    3.  Calculates `s_k = tempPhi - phi` (change in `phi`).
    4.  Calculates `y_k = tempdLdPhi - dLdPhi` (change in gradient `dL/dPhi`).
    5.  Stores `s_k` and `y_k` in the circular buffers `BFGSs` and `BFGSy` at slot `BFGSPT`.
    6.  Updates `BFGSPT` to point to the next slot.
  - **Arguments:**
    - `phi :: REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:)`: `phi` before the line search.
    - `tempPhi :: REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:)`: `phi` after the line search.
    - `dEdPhi :: REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:)`: `dE/dphi` at `tempPhi`.
    - `dLdPhi :: REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:)`: `dL/dphi` at `phi` (before line search).
    - `totEleNum(:) :: REAL(KIND=DP), INTENT(IN)`: Total number of electrons per spin channel.
    - `dV :: REAL(KIND=DP), INTENT(IN)`: Real-space grid volume element.

- **`SUBROUTINE ChemicalPotential(dEdPhi, phi, totEle, mu0)`**:
  - **Description:** Calculates the chemical potential `mu0` for each spin channel. The formula used is `mu0 = integral( (dE/dPhi_spin) * Phi_spin ) / (2 * N_electrons_spin)`. This ensures that the gradient `dL/dPhi = dE/dPhi - 2*mu*Phi` respects the electron number constraint when `d L / d mu = 0`. It performs an MPI sum reduction for `mu0`.
  - **Arguments:**
    - `dEdPhi :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: `dE/dPhi`.
    - `phi :: REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN)`: Current `phi`.
    - `totEle(:) :: REAL(KIND=DP), INTENT(IN)`: Total electrons per spin.
    - `mu0(:) :: REAL(KIND=DP), INTENT(INOUT)`: Calculated chemical potential(s).

# Important Variables/Constants

- **`MBFGS`**: Defines the memory/history length for the L-BFGS algorithm.
- **`BFGSs`, `BFGSy`, `BFGSrho`**: Core arrays storing the history required for the L-BFGS two-loop recursion.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP`.
- **`OutputFiles`**: Module is used, but `outputUnit` is not directly referenced in these subroutines.
- **`MPI_Functions`**: For `ReduceRealLevel1` (used in `ChemicalPotential` and via `InnerProd`).
- **`RhoDirCG`**: For the `InnerProd` function, which is used in `BFGSDirection` to compute scalar products like `s_i^T q` and `y_i^T r`.
- **`CellInfo`**: (Implicit via `InnerProd` from `RhoDirCG`) For `cell%dV`.

This module provides the logic for the L-BFGS search direction and history update. It would be called by a higher-level electron density optimization routine (e.g., `RhoOptN` or a dedicated BFGS optimizer module). That optimizer would handle the overall iteration loop, perform line searches along the direction provided by `BFGSDirection`, and then call `UpdateBFGS` with the results of the line search. The `ChemicalPotential` subroutine is crucial for ensuring electron number conservation in the optimization.
