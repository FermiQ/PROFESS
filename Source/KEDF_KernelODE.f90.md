# Overview

The `IntKernelODE` module is responsible for numerically solving the Ordinary Differential Equations (ODEs) that define the non-local kernels used in various Kinetic Energy Density Functionals (KEDFs). Specifically, it is used to generate the kernel function `w(eta)` and its derivatives, which are essential components for KEDFs like Wang-Govind-Carter (WGC) and Corrected Average Thomas-Fermi (CAT).

The module interfaces with an external Runge-Kutta solver package (referred to as `rksuite.f` in comments) to perform the numerical integration of these ODEs. The ODEs are typically solved by integrating from `eta = infinity` (where an asymptotic value `wInf` is known) down to `eta = 0`. The results of the integration – `w(eta)` and its derivatives – are stored in module-level arrays for subsequent use by the respective KEDF modules (e.g., `KEDF_WGC`, `KEDF_CAT`, `KEDF_HC10`).

# Key Components

- **`MODULE IntKernelODE`**: The main container for ODE solving routines.

- **`SUBROUTINE makeKernel_ode1`**:
  - **Description:** Solves a **1st order** ODE for a kinetic energy kernel. It integrates the ODE defined in the `SingleDep` subroutine to obtain `w(eta)` (stored in `w`) and `w'(eta)` (stored in `wp`).
  - **Method:** Uses the external `rksuite` solver, integrating from a large `eta` down to `eta=0`.
  - **Output:** Populates the module-level arrays `eta`, `w`, and `wp`.

- **`SUBROUTINE makeKernel2`**:
  - **Description:** Solves a **2nd order** ODE for a kinetic energy kernel. It can solve for different KEDF types based on the `KEDFtype` flag (1 for CAT, 2 for WGC). It integrates the ODE defined in either the `CAT` or `WGC` subroutine.
  - **Method:** Uses the external `rksuite` solver, integrating from a large `eta` down to `eta=0`.
  - **Output:** Populates the module-level arrays `eta`, `w` (for `w(eta)`), `w1` (for `w'(eta)`), and `w2` (for `w''(eta)`).

- **`SUBROUTINE WGC(T, Y, YP)`**:
  - **Description:** Defines the 2nd order ODE for the Wang-Govind-Carter (WGC) KEDF kernel. This subroutine is passed to the `rksuite` solver.
  - **Arguments (for solver):**
    - `T :: REAL(KIND=DP)`: Input, current value of `eta`.
    - `Y :: REAL(KIND=DP), DIMENSION(2)`: Input, `Y(1) = w(eta)`, `Y(2) = w'(eta)`.
    - `YP :: REAL(KIND=DP), DIMENSION(2)`: Output, `YP(1) = w'(eta)`, `YP(2) = w''(eta)`.

- **`SUBROUTINE CAT(T, Y, YP)`**:
  - **Description:** Defines the 2nd order ODE for the Corrected Average Thomas-Fermi (CAT) KEDF kernel. This subroutine is passed to the `rksuite` solver.
  - **Arguments:** Same structure as `WGC`.

- **`SUBROUTINE SingleDep(T, Y, YP)`**:
  - **Description:** Defines a 1st order ODE for a "single density dependent kernel" (e.g., used in HC10 KEDF). This subroutine is passed to the `rksuite` solver (via `makeKernel_ode1`).
  - **Arguments (for solver):**
    - `T :: REAL(KIND=DP)`: Input, current value of `eta`.
    - `Y :: REAL(KIND=DP), DIMENSION(1)`: Input, `Y(1) = w(eta)`.
    - `YP :: REAL(KIND=DP), DIMENSION(1)`: Output, `YP(1) = w'(eta)`.

- **`SUBROUTINE Clean`**:
  - **Description:** Deallocates the module-level arrays (`eta`, `w`, `w1`, `w2`, `wp`) that store the kernel solution.

# Important Variables/Constants

- **Module-Level Parameters for ODEs:**
    - `ode_alpha :: REAL(KIND=DP)`: Parameter `alpha` for the kernel ODE (e.g., used in CAT).
    - `ode_beta :: REAL(KIND=DP)`: Parameter `beta` for the kernel ODE.
    - `ode_gamma :: REAL(KIND=DP)`: Parameter `gamma` for the kernel ODE (e.g., used in CAT).
    - `wInf :: REAL(KIND=DP)`: The asymptotic value of `w(eta)` as `eta -> infinity`, used as the initial condition for the ODE integration.
    - `KEDFtype :: INTEGER`: An integer flag that selects which kernel ODE to solve in `makeKernel2` (1 for `CAT`, 2 for `WGC`).

- **Module-Level Arrays for Kernel Storage:**
    - `eta(:) :: REAL(KIND=DP), ALLOCATABLE`: Stores the grid of `eta` points at which the solution is evaluated.
    - `w(:) :: REAL(KIND=DP), ALLOCATABLE`: Stores the kernel function `w(eta)`.
    - `w1(:) :: REAL(KIND=DP), ALLOCATABLE`: Stores the first derivative `w'(eta)` (typically from 2nd order ODEs).
    - `w2(:) :: REAL(KIND=DP), ALLOCATABLE`: Stores the second derivative `w''(eta)` (typically from 2nd order ODEs).
    - `wp(:) :: REAL(KIND=DP), ALLOCATABLE`: Stores the first derivative `w'(eta)` (typically from 1st order ODEs via `makeKernel_ode1`).

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`Constants`**: For `DP` (double precision kind parameter).
- **`Output`**: For `Wrtout` (writing messages).
- **External ODE Solver (`rksuite.f`)**: This is a major external dependency. The module uses subroutines from this package, such as:
    - `SETUP_RK`: To initialize the solver.
    - `UT_RK`: To perform an integration step and get the solution at a desired point.
    - `RK_STATISTICS` (referenced as `STAT` in `EXTERNAL` statement): To get statistics about the integration process.
  The actual ODEs to be solved (`WGC`, `CAT`, `SingleDep`) are passed as arguments to `UT_RK`.

**Workflow:**
A KEDF module that requires a non-local kernel (e.g., `KEDF_WGC`, `KEDF_CAT`, `KEDF_HC10`) will typically:
1. Set the appropriate parameters in `IntKernelODE` (`ode_alpha`, `ode_beta`, `ode_gamma`, `wInf`, `KEDFtype`).
2. Call either `makeKernel_ode1` or `makeKernel2`.
3. This called routine then uses the external `rksuite` solver, providing one of `WGC`, `CAT`, or `SingleDep` as the function defining the ODE system.
4. The solver populates the `eta`, `w`, `w1`, `w2`, or `wp` arrays.
5. The calling KEDF module then uses these arrays (often after splining them onto its own G-space grid, as seen in `KEDF_HC10`'s `FillCAT` or `KEDF_WGC`'s `FillWGC`) to construct its kinetic energy and potential.
6. At the end of the program, `Clean` is called to deallocate the stored kernel arrays.
