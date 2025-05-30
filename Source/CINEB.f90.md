# Overview

This is a Climbing Image Nudged Elastic Band (CINEB) code that determines transition states using an Orbital Free DFT (OFDFT) code. It manages ion movements and performs single point energy/force calculations by calling an external OFDFT program.

# Key Components

The program `CINEB` orchestrates the entire calculation. Key subroutines and functions include:

**Input Subroutines:**
- `SUBROUTINE ReadCinebInputFile`: Reads the main input file (`CINEB.inp`) containing cell parameters, atom specifications, image coordinates, and frozen atom constraints.
- `SUBROUTINE ReadCinebParamFile`: Reads the parameter file (`CINEB.param`) for settings like `nudgeType`, `initialConditions`, and `forceCut`.

**Calculation Subroutines:**
- `SUBROUTINE TransitionStateVelocityVerlot()`: Implements a velocity Verlet-like algorithm (described as QUICKMIN in comments) to minimize the total force on the band of images, searching for the minimum energy path.
- `SUBROUTINE TotalCalculation(numCalcs, cell, cartPosition, atomicSymbol, atomicPSP, trueForce, springForce, totalForce, energy, stressTensor, springConstant)`: Calculates all relevant forces (true DFT forces, spring forces between images, and the total combined force) and energies for the current set of image positions.
- `SUBROUTINE TrueCalculation(numcalcs, cell, cartPosition, atomicSymbol, atomicPSP, energy, trueForce, stressTensor)`: Manages the execution of the external OFDFT program for each image to get energies, true forces, and stress tensors. It prepares input files for OFDFT, runs the calculations (potentially in parallel directories), and harvests the results.
- `SUBROUTINE GetSpringConstant(energy, trueForce, massPosition, springConstant)`: Determines appropriate spring constants between images, potentially varying them based on image energies, to ensure a smooth band.
- `SUBROUTINE GetTotalForce(positions, trueForces, springForces, energies, springConstant, totalForces)`: Calculates the total effective force on each ion in each image, by combining the true forces from DFT and the spring forces, according to the selected `nudgeType` (NEB, CINEB, etc.). It also handles the climbing image logic.
- `SUBROUTINE HenkelTangent(positions, energies, tangents)`: Calculates tangents for the elastic band using the improved tangent estimate method by Henkelman et al. (J. Chem. Phys. 113, 9978 (2000)).
- `SUBROUTINE HenkelForce(positions, tangents, trueForces, springForces, energies, springConstant, totalForces)`: Calculates the spring and total forces based on the Henkelman tangent definition. Implements the climbing image force modification if `nudgeType` corresponds to CINEB with this tangent.
- `SUBROUTINE SplineTangent(positions, tangents)`: Calculates tangents by fitting a cubic spline to the path defined by the images.
- `SUBROUTINE SplineForce(tangents, trueForces, springForces, energies, totalForces)`: Calculates the spring and total forces based on the spline-derived tangents. Implements the climbing image force modification if `nudgeType` corresponds to CINEB with this tangent.

**Mathematical Subroutines:**
- `SUBROUTINE Spline(xvalues, yvalues, spfit)`: Performs a cubic spline interpolation. (Older version)
- `SUBROUTINE Spline2(x, y, spfit)`: Gives the coefficients for a cubic spline fit going through data points (x,y) with natural boundary conditions. Solves using Crout factorization. (Newer version by Linda Hung)
- `SUBROUTINE Recenter(cell, cartPosition, atomicSymbol, massPosition)`: Calculates ion positions in center-of-mass coordinates to prevent net translation of the images. It also includes logic for rotational alignment (though the `Rotate` call is commented out in the version provided).
- `SUBROUTINE Rotate(vector, theta)`: Rotates a 3D vector by specified angles around x, y, and z axes.
- `FUNCTION CartToFract(cell, cartPosition)`: Converts Cartesian coordinates to fractional coordinates based on the cell lattice vectors.
- `FUNCTION FractToCart(cell, fractPosition)`: Converts fractional coordinates to Cartesian coordinates.
- `SUBROUTINE GaussianElimination(NR, NCLHS, LEFT, NCRHS, RIGHT)`: Performs Gaussian elimination with partial pivoting to solve systems of linear equations.

**Output Subroutines:**
- `SUBROUTINE OutputWriteHeader(cell)`: Writes the header section to the main output file (`cineb.out`), including cell parameters and column titles for the iteration summary.
- `SUBROUTINE OutputWrite(cartPosition, massPosition, energy, trueForce, springForce, totalForce, stressTensor, numcalcs, numiterations, linecounts, badStep, timeStep)`: Writes detailed information for the current iteration to `cineb.out`, including energies, forces for each image, and convergence metrics.
- `SUBROUTINE OutputWriteFooter()`: Writes a footer to `cineb.out` at the end of the calculation.
- `SUBROUTINE WriteRestartInformation(cartPosition, cell, atomicPSP, atomicSymbol)`: Writes current image positions and other relevant data to `restart.out`, allowing the calculation to be resumed.

# Important Variables/Constants

- `numAtoms :: INTEGER`: The number of atoms in the simulation cell.
- `numAtomTypes :: INTEGER`: The number of different atomic species.
- `numImages :: INTEGER`: The number of images (intermediate states) used to define the path between initial and final states in the NEB calculation.
- `nudgeType :: INTEGER`: Controls the NEB algorithm:
    - `-1`: No band (likely for single point calculations on endpoints).
    - `0`: Plain Elastic Band.
    - `1`: NEB with "improved" tangent (Henkelman).
    - `2`: CINEB with "improved" tangent (Henkelman), allows climbing image.
    - `3`: NEB with original spline tangent.
    - `4`: CINEB with original spline tangent, allows climbing image.
    - `5`: NEB with new spline tangent (Spline2).
    - `6`: CINEB with new spline tangent (Spline2), allows climbing image.
- `initialConditions :: INTEGER`: Defines how initial image configurations are handled:
    - `1`: Use provided initial coordinates directly.
    - `2`: Relax the first and last images using an optimization method (QUICKMIN mentioned) before starting the NEB.
- `forceCut :: REAL(KIND=DP)`: Convergence criterion for the maximum force on any atom in the band. The NEB calculation stops when the maximum force falls below this value.
- `cell :: REAL(KIND=DP), DIMENSION(3,3)`: Lattice vectors defining the simulation cell.
- `cartPosition :: REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE`: Stores the Cartesian coordinates of atoms for each image (`cartPosition(image, atom_index, coordinate_index)`).
- `atomicSymbol :: CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE`: Stores the chemical symbols for each atom (e.g., "Si", "Al").
- `atomicPSP :: CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE`: Stores the filenames or identifiers for the pseudopotentials of each atom type.
- `frozenAtoms :: INTEGER, DIMENSION(:,:), ALLOCATABLE`: A 2D array (`frozenAtoms(atom_index, coordinate_index)`) where a value of 0 indicates the atom's movement in that Cartesian direction (x,y,z) is constrained (frozen).
- `posType :: CHARACTER(LEN=1)`: Indicates the format of input atomic positions: 'F' for fractional coordinates, 'C' for Cartesian.
- `transOutputUnit, transImagesOutputUnit, transImagesOutputForcesUnit :: INTEGER`: Fortran unit numbers for various output files.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

**Internal Module Dependencies:**
- `USE CONSTANTS, ONLY : DP, pi, AtomicMass`: Uses double precision kind (`DP`), the constant `pi`, and `AtomicMass` table from the `CONSTANTS` module.
- `USE MATHFUNCTIONS, ONLY : Inverse, Volume, Cross, Vecmul`: Uses matrix inverse (`Inverse`), cell volume calculation (`Volume`), vector cross product (`Cross`), and matrix-vector multiplication (`Vecmul`) from the `MATHFUNCTIONS` module.

**External Program Interactions:**
- The program calls an external OFDFT executable named `PROFESS` using the `SYSTEM` intrinsic. It constructs input files for `PROFESS` (e.g., `IMAGE*/IMAGE*.inpt`, `IMAGE*/IMAGE*.ion`) in separate directories for each image and then parses output files (e.g., `IMAGE*/IMAGE*.trans`, `IMAGE*/IMAGE*.force.out`) to retrieve energies and forces.

**File Interactions:**
- **Input files:**
    - `CINEB.inp`: Main input file specifying geometry, images, and atom types.
    - `CINEB.param`: Contains parameters controlling the NEB calculation.
    - `ORBITALFREE.inp`: Template input file for the OFDFT calculations. Copied and modified for each image.
    - `*.recpot`, `*.den` (potentially): Pseudopotential and density files, copied into image directories for OFDFT calculations.
- **Output files:**
    - `cineb.out`: Main human-readable output summarizing the NEB progression.
    - `IMAGE*/IMAGE*.log`: Log file from the `PROFESS` OFDFT calculation for each image.
    - `IMAGE*/IMAGE*.trans`: File from OFDFT containing energy, forces, and stress for an image.
    - `IMAGE*/IMAGE*.force.out`: Another output from OFDFT, its specific content related to forces.
    - `MONITOR.dat`: Used to track the status of parallel OFDFT jobs.
    - `executeOFDFT.s`: A shell script generated on the fly to run OFDFT calculations for all images.
    - `restart.out`: Contains necessary information (cell, atomic positions for all images) to restart the CINEB calculation.
- **Temporary/Working files:**
    - Files within `IMAGE*` directories (e.g., `IMAGE*.inpt`, `IMAGE*.ion`, `IMAGE*.calc.ion`).

The program manages a series of calculations across multiple "images" forming a path. It iteratively refines these images to find a minimum energy path and identify transition states. The core loop involves:
1.  Calculating true forces via external OFDFT calls (`TrueCalculation`).
2.  Determining inter-image spring forces.
3.  Calculating tangents to the path.
4.  Projecting out components of forces and applying nudging/climbing image logic (`GetTotalForce`).
5.  Updating image positions based on these forces (`TransitionStateVelocityVerlot`).
6.  Writing output and checking for convergence.
