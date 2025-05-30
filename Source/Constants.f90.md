# Overview

The `Constants` module serves as a centralized definition point for a wide range of universal physical constants, conversion factors, and mathematical constants. These values are intended to be static. The module also includes a utility function, `AtomicMass`, to retrieve the atomic mass for a given chemical element symbol. The selection of constants aims for consistency with those used in other established codes like CASTEP, though periodic checks against updated standard values might be necessary.

# Key Components

- **`MODULE Constants`**:
  - The main container for all defined constants and the `AtomicMass` function.

- **`FUNCTION AtomicMass(symbol) RESULT(mass_value)`**:
  - **Description:** Takes a chemical element symbol (e.g., "H", "Si", "Al") as input and returns its atomic mass in atomic units (where the mass of a carbon-12 atom is 12). The function uses a `SELECT CASE` structure to look up the mass based on the trimmed input symbol. If the symbol is not found, it stops execution with an error message.
  - **Arguments:**
    - `symbol :: CHARACTER(LEN=3), INTENT(IN)`: The chemical symbol of the element. The function internally trims whitespace.
  - **Return Value:**
    - `AtomicMass :: REAL(KIND=DP)`: The atomic mass of the specified element in atomic units.

# Important Variables/Constants

All constants are defined as `PARAMETER`s, meaning they are fixed at compile time.

**Precision and Length Parameters:**
- `DP :: INTEGER, PARAMETER`: Defines the kind parameter for double precision floating-point numbers (set to 8). This ensures approximately 14 significant digits.
- `systemNameLen :: INTEGER, PARAMETER`: Maximum length for system names or similar character strings (set to 80).

**Mathematical and General Constants:**
- `tiny :: REAL(kind=DP), PARAMETER`: A very small number (1.0E-12_DP).
- `pi :: REAL(kind=DP), PARAMETER`: The mathematical constant pi (3.14159...).
- `golden :: REAL(kind=DP), PARAMETER`: The golden ratio (0.38196...).
- `imag :: COMPLEX(KIND=DP), PARAMETER`: The imaginary unit, i (0.0 + 1.0i).

**Physical Constants:**
- `qElec :: REAL(kind=DP), PARAMETER`: Elementary charge (1.602176462E-19 C).
- `light :: REAL(kind=DP), PARAMETER`: Speed of light in vacuum (299792458.0 m/s).
- `bohr :: REAL(kind=DP), PARAMETER`: Bohr radius (0.529177208607388 Angstroms).
- `emass :: REAL(kind=DP), PARAMETER`: Electron mass (9.1093826E-28 g).
- `mol :: REAL(kind=DP), PARAMETER`: Avogadro's number (6.02214179E23 mol^-1).
- `boltzmann :: REAL(kind=DP), PARAMETER`: Boltzmann constant (3.16682E-6 Hartree/K).
- `mu0 :: REAL(kind=DP), PARAMETER`: Permeability of vacuum (4*pi*1.0E-7 N/A^2).
- `eps0 :: REAL(kind=DP), PARAMETER`: Permittivity of vacuum (derived from `mu0` and `light`).

**Conversion Factors:**
- `mToA :: REAL(kind=DP), PARAMETER`: Meters to Angstroms (1.0E10).
- `rydbergToeV :: REAL(kind=DP), PARAMETER`: Rydberg to electron-volts (13.60569193 eV).
- `hartreeToeV :: REAL(kind=DP), PARAMETER`: Hartree to electron-volts (derived, approx 27.211 eV).
- `fundamentalTime :: REAL(kind=DP), PARAMETER`: Seconds per atomic unit of time (2.41888...E-17 s/a.u.).
- `auPressure :: REAL(kind=DP), PARAMETER`: Atomic units of pressure to Pascals (1.0 / 2.9421912E13 a.u./Pa).
- `auToGPa :: REAL(kind=DP), PARAMETER`: Atomic units of pressure to Gigapascals.

**Other Variables:**
- `machPrec :: REAL(KIND=DP)`: Declared to hold machine precision. *Note: This variable is declared but not initialized within the module as a `PARAMETER`. Its value would need to be set at runtime or be a `PARAMETER` itself if intended as a fixed constant.*

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

The `Constants` module is designed to be a foundational module with no dependencies on other custom modules within the PROFESS codebase for its definitions. It is, however, extensively used by virtually all other modules that require physical constants, conversion factors, or atomic masses.
