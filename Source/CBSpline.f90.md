# Overview

Used in spline method for ion-ion and ion-electron interactions

# Key Components

- `SUBROUTINE GetCardinalBSpline(spline, order, x)`: Solve for M_n(x) for integer values of x, 0<x<n, where M_n is the nth order Cardinal B-spline.  Use iterative method with M_n(x) = x/(n-1)*M_(n-1)(x) + (n-x)/(n-1)*M_(n-1)(x-1)
- `SUBROUTINE FillB(b1, b2, b3, order, totX, tot3, locOffset)`: This routine is during module setup.  It is called only once, assuming fixed recip grid size/spacing and b-spline order
- `SUBROUTINE BSplineProduct(array, results, squareProduct)`: Muliplies an array with b-splines, which were previously set up.

# Important Variables/Constants

- `bSpline1, bSpline2, bSpline3`: Allocatable complex arrays used to store B-spline values for each dimension.
- `splineOrder :: INTEGER`: Defines the order of cardinal b-splines, typically used for ion-ion and ion-electron interaction terms. Default value is 10.

# Usage Examples

[If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.]

# Dependencies and Interactions

This module uses entities from the following modules:
- `Constants`: Specifically, `DP` (double precision kind parameter), `PI` (mathematical constant pi), and `IMAG` (imaginary unit).
- `OutputFiles`: Used for error handling, specifically the `ERROR` subroutine is called.

The module provides the following key subroutines for use by other parts of the system: `GetCardinalBSpline`, `FillB`, and `BSplineProduct`. It also defines module-level allocatable arrays `bSpline1`, `bSpline2`, `bSpline3` and an integer `splineOrder`.
