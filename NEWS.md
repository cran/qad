
<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# qad 1.0.1

## Changes

-   Added a vignette
-   Added automatic testing procedures

# qad 1.0.0

## Main changes

-   Rcpp code implementation for the main function qad()
-   Code improvements for zeta1() and ECBC() which replaces the
    deprecated function emp_c\_copula()
-   Adapted permutation test for independence which works irrespectively
    of the marginal distribution
-   Added a new (experimental) bootstrap test which tests for asymmetry
    in dependence (thus provides a p-value for the measure a)

## Minor changes

-   Additional parameters for graphics (title, x.axis, y.axis)
-   Added warnings for checkerboard resolution less or equal 3
-   the deprecated function cci() was finally removed

# qad 0.2.0

## Main changes

-   Implementation of a new p-value (based on Monte Carlo simulation).
-   New implementation of the distribution function and quantile
    function of qad in the setting of independence: pqad(), qqad().
-   Function cci() is deprecated.
-   Additional parameter configurations in pairwise.qad: (1) double 0
    entries can be removed, (2) a minimum necessary resolution can be
    selected.

## Bug Fixes

-   Correction of a wrong p-value (by permutation test) of qad for a
    checkerboard resolution of N = 1.
-   Code simplifications (replace several functions of the data.table
    package by dplyr functions).

# qad 0.1.2

## Internals

-   Update of functions using melt from data.table

# qad 0.1.1

## Bugfixes and changes

-   New color options are available in the plot functions
-   New labeling in the heatmap plot
-   Updated output of the main function qad
-   New option for parallel computing

## Internals

-   Package news are now provided as markdown document (NEWS.md)
-   Package description in form of a README.md file is available
