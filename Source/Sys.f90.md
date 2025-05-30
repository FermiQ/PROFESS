# Overview

The `SYS` module (defined in the file `System.f90`) acts as a central repository for many of the global state variables and core data arrays that describe the physical system being simulated. This includes the primary quantities like electron density (`rhoR`), real-space potential (`potReal`), ionic forces (`forceIon`), the system stress tensor (`stress`), and a breakdown of energy components (`energy`).

It also holds various simulation parameters and flags, such as information about frozen ions (`frozenIon`), reference densities for KEDFs (`rho0`, `rhoS`), and flags controlling their behavior (`bvac`, `hold0`, `holdS`). The module provides a main setup routine (`SetupSystem`) to initialize these arrays and parameters based on inputs (like energy cutoff and KEDF choice) and a cleanup routine (`CleanSystem`) to deallocate them.

The comment "This information is limited to be known only by a SELECT FEW modules" suggests an attempt to control global variable access, primarily to `Initializer`, `Optimizer`, and `Output` modules.

# Key Components

- **`MODULE SYS`**: The main container for global system variables and setup routines.

- **`TYPE gridPack`**:
  - **Description:** A derived data type intended for multi-grid approaches (comment mentions "AMD" - Adaptive Mesh Refinement).
  - **Fields:**
    - `rhoCore :: REAL(KIND=DP), DIMENSION(:,:,:), POINTER`: Pointer to a core density on a grid.
    - `section :: REAL(KIND=DP), DIMENSION(3,3)`: Describes fractional sections of the grid for multi-grid schemes.

- **`SUBROUTINE SetupSystem(energyCutoff, kinetic)`**:
  - **Description:** Initializes the core data arrays and parameters for the system.
    1.  Calculates machine precision.
    2.  Calls `InitializeFFT` (from `SetupFFT`) to determine FFT grid dimensions (`m1G`, `m2G`, `m3G` globally; `n1G`, `n2G`, `n3G`, `n3Goff` locally) based on `energyCutoff` or `gridSpacing`.
    3.  Sets up MPI information for ion distribution if running in parallel (`numIonLoc`, `numIonInit`).
    4.  Defines real-space grid boundaries (`xmin` to `xmax`, etc.).
    5.  Allocates main system arrays: `rhoR`, `interior`, `potReal`, `forceIon`.
    6.  Initializes `rhoR` to a uniform density based on the total number of electrons derived from `cell%elementTable%chargeTot`. Handles spin polarization by adjusting density for `magmom`.
    7.  Allocates and initializes the `grids` array (of `gridPack` type).
  - **Arguments:**
    - `energyCutoff :: REAL(KIND=DP), INTENT(IN)`: Kinetic energy cutoff for determining grid size.
    - `kinetic :: INTEGER, INTENT(IN)`: Parameter indicating the choice of KEDF (influences `kePotDim`, though `kePotDim` is not used further in this routine).

- **`SUBROUTINE CleanSystem(ionKeepOpt)`**:
  - **Description:** Deallocates the major system arrays: `rhoR`, `interior`, `potReal`, `forceIon`, `grids`. Optionally, if `ionKeepOpt` is `.FALSE.` (default), it also deallocates pseudopotential data stored within `cell%elementTable`. Deallocates `frozenIon` if it was allocated. Calls `CleanFFT` (from `Fourier` module).
  - **Arguments:**
    - `ionKeepOpt :: LOGICAL, INTENT(IN), OPTIONAL`: If `.TRUE.`, pseudopotential data within `cell` is not deallocated.

- **`FUNCTION GetMachinePrecision() RESULT(eps)`**:
  - **Description:** Calculates and returns the machine precision for `REAL(KIND=DP)` numbers.

- **`SUBROUTINE CalculateElectronNumber(cellLoc, rho, dimX, dimY, dimZ, dimTotal, totEleNum)`**:
  - **Description:** Calculates the total number of electrons by integrating a given density `rho` over the volume of `cellLoc`. `dimTotal` is the total number of global grid points for calculating `dV`. Performs an MPI sum reduction for `totEleNum`.
  - **Arguments:** `cellLoc` (cell data), `rho` (density array), `dimX,Y,Z` (dims of `rho`), `dimTotal` (global grid points), `totEleNum` (OUT).

- **`SUBROUTINE CleanSystemForESPDOpt()`**:
  - **Description:** A variant of `CleanSystem`, seemingly for a specific "ESPD" optimization. It deallocates `rhoR`, `interior`, and `frozenIon`, and mentions other variables (`rhocore`, `auxidm`, `fixingvariables`) that should be deallocated but are not standard module variables, indicating this might be an incomplete or specialized cleanup. Calls `CleanFFT`.

# Important Variables/Constants (Module-Level Global Data)

- **Energy & Magnetism:**
    - `energy(9) :: REAL(KIND=DP)`: Array storing the total energy (`energy(1)`) and its various components (kinetic, external, Coulomb, XC, ion-ion, TF, vW, Nonlocal KEDF). Initialized to `0.0_DP`.
    - `magmom :: REAL(KIND=DP)`: Total magnetic moment of the system, used for initializing spin-polarized densities (default `0.0_DP`).
- **Core Electronic Structure Arrays (Allocatable):**
    - `rhoR(:,:,:,:) :: REAL(KIND=DP), TARGET`: The real-space electron density. Dimensions are `(n1G, n2G, n3G, numSpin)`.
    - `interior(:,:,:,:) :: LOGICAL`: A mask, likely indicating grid points within the primary cell versus boundary/ghost cells (though initialized to all `.TRUE.` for periodic conditions).
    - `potReal(:,:,:,:) :: REAL(KIND=DP)`: The total effective real-space potential.
- **Ionic & Cell Variables (Allocatable or direct):**
    - `frozenIon(:,:) :: LOGICAL`: Mask indicating which ionic degrees of freedom (atom index, Cartesian direction) are constrained.
    - `forceIon(:,:,:) :: REAL(KIND=DP)`: Forces on the ions. Dimensions `(numIon, 3, numForceComponents)`. `forceIon(:, :, 1)` is total force.
    - `stress(3,3) :: REAL(KIND=DP)`: The system stress tensor.
- **Multi-Grid & KEDF Parameters:**
    - `grids(:) :: TYPE(gridPack)`: Array for storing `gridPack` types, for multi-grid methods.
    - `gridSpacing :: REAL(KIND=DP)`: If `energyCutoff` is not used, this defines grid density (default -1.0).
    - `bvac :: LOGICAL`: Flag for vacuum corrections in some KEDFs (default `.FALSE.`).
    - `rho0 :: REAL(KIND=DP)`: Reference density for WT/WGC KEDFs (default -1.0).
    - `rhoS :: REAL(KIND=DP)`: Reference density for WGC kernel expansion (default -1.0).
    - `hold0, holdS :: LOGICAL`: Flags indicating if `rho0` and `rhoS` are fixed during cell optimization (affect stress calculation). Defaults `.FALSE.`.
    - `LumgpExp, LumgpFactor :: REAL(KIND=DP)`: Parameters for the MGP KEDF (defaults 0.0).
- **Parallelism Specific (Compiled if `__USE_PARALLEL` is defined):**
    - `numIonLoc :: INTEGER`: Number of ions handled by the local processor.
    - `numIonInit :: INTEGER`: Starting global index of ions handled by the local processor.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- **`CONSTANTS`**: For `DP` and `machprec` (though `machprec` is calculated locally too).
- **`MPI_Functions`**: For MPI-related variables (`rankGFFT`, `sizeGFFT`) and utilities (`ReduceRealLevel1`, `TITLE`).
- **`OutputFiles`**: For `outputUnit`.
- **`CellInfo`**: For grid dimensions (`m1G` etc.), `cell` derived type, and `numSpin`.
- **`SetupFFT`**: Critically for `InitializeFFT` which determines grid dimensions and sets up FFT plans.
- **`Fourier`**: The `CleanSystem` routine calls `CleanFFT` (from the `Fourier` shim module, which would call the `_NEW` version).

The `SYS` module is central to PROFESS, holding most of the key physical quantities and simulation parameters that are shared across different computational stages. Its `SetupSystem` routine is a major part of the initialization process.
