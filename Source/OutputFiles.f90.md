# Overview

The `OutputFiles` module serves two primary purposes:
1.  It defines global `PARAMETER` constants for standard Fortran unit numbers to be used for general program output (`outputUnit`) and error messages (`errorUnit`).
2.  It provides a standardized way to report errors and terminate the program through the generic `Error` interface (`OneError` for single messages, `ManyError` for multiple messages).

Error reporting is designed to be MPI-aware: only MPI rank 0 (as determined by the `outRank` variable) will print the error messages before the program stops. This prevents redundant error messages from all processes in a parallel run.

# Key Components

- **`MODULE OutputFiles`**: The main container module.

- **Module-Level Parameters:**
    - `outputUnit :: INTEGER, PARAMETER`: Defines the Fortran unit number for standard program output (default: 8).
    - `errorUnit :: INTEGER, PARAMETER`: Defines the Fortran unit number for error messages (default: 9). This unit is used by `OneError` to also write the error message, ensuring it goes to the designated error stream.

- **Module-Level Variables:**
    - `message :: CHARACTER(LEN=500)`: A character string declared at the module level. While present, it is not directly used by the `Error` interface subroutines, as they take their messages as arguments. It might be intended for other general-purpose messaging not shown or for future use.
    - `outRank :: INTEGER`: Stores the MPI rank of the current process. This variable is expected to be set by an MPI initialization routine (e.g., in the `MPI_Functions` module) to `rankGlobal`. Its default value is 0, making error reporting work correctly in serial runs or if MPI is not used/initialized.

- **`INTERFACE Error`**:
  - **Description:** A generic interface allowing a call to `CALL Error(...)` to be resolved to either `OneError` or `ManyError` based on the arguments provided.
  - **Module Procedures:**
    - `OneError(unitOut, msg)`
    - `ManyError(unitOut, msg, num)`

- **`SUBROUTINE OneError(unitOut, msg)`**:
  - **Description:** Handles the reporting of a single error message. If `outRank` is 0, it prints a formatted error block containing the string `msg` to the specified `unitOut` and also to the globally defined `errorUnit`. After printing, it calls `STOP` to terminate program execution.
  - **Arguments:**
    - `unitOut :: INTEGER, INTENT(IN)`: The Fortran unit number to which the primary error message is written (e.g., standard output).
    - `msg :: CHARACTER(LEN=*), INTENT(IN)`: The error message string.

- **`SUBROUTINE ManyError(unitOut, msg, num)`**:
  - **Description:** Handles the reporting of multiple error messages stored in an array. If `outRank` is 0, it prints a formatted error block. It then iterates through the `msg` array (up to `num` messages, or 1 if `num` is not provided) and prints each message to `unitOut`. After printing, it calls `STOP`.
  - **Arguments:**
    - `unitOut :: INTEGER, INTENT(IN)`: The Fortran unit number for output.
    - `msg(:) :: CHARACTER(LEN=*), INTENT(IN)`: An array of error message strings.
    - `num :: INTEGER, OPTIONAL, INTENT(IN)`: The number of messages from the `msg` array to print. Defaults to 1 if not present.

# Important Variables/Constants

- **`outputUnit`, `errorUnit`**: Define the standard channels for program output and error logging.
- **`outRank`**: Crucial for controlling output in parallel MPI runs, ensuring that messages (especially error messages leading to termination) are printed only once by the master process.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

- This module itself has no `USE` statements pointing to other custom modules from the PROFESS project. It is a foundational utility module.
- **Implicit Dependency on MPI Setup:** For parallel runs, the correct functioning of the rank-dependent output in the `Error` subroutines relies on `outRank` being properly set by an MPI initialization routine (e.g., `MPI_Functions.InitializeMPI()`).
- Other modules throughout the application will `USE OutputFiles` to access the standard `outputUnit` and `errorUnit` for their I/O operations, and to call the `Error` interface for standardized program termination upon encountering fatal errors.

The `Error` interface provides a consistent way to handle fatal errors, ensuring that a descriptive message is logged before the program halts.
