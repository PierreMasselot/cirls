# Current development version

## New features
- New method `shapeConstr` to build constraint matrices for shape-constrained splines. Currently works with classes `ns`, `bs` (from package `splines`), `ps` and `onebasis` (from `dlnm`). Also includes a default method for more general bases.
- New function `edf` to compute observed and expected degrees of freedom for a fitted `cirls` object.
- Method `logLik.cirls` for AIC and BIC computation.
- New function `uncons` to return the unconstrained model.

## Changes
- Changed the default QP solver to `quadprog` after some expriments.
- Added the argument `complete` in inferential functions, to allow keeping or discarding aliased coefficients. Same interpretation as in `vcov.lm`.
- `vcov` now allows returning the usual variance-covariance returned by `vcov.glm` when `constrained = FALSE`.
- `simu_coef` renamed as `simulCoef` and now includes an argument to set the seed. 
- `check_cmat` has been renamed `checkCmat`.
- Now the `Cmat`, `lb` and `ub` used are not returned in the `control` object from the result of `glm` with `cirls.fit`.
- Now `lb` and `ub` can be passed by term.
- The element `aic` of a fitted `cirls` object is penalised by the number of active constraints.

## Bug fixes
- Fixed issue with R matrix when there was less observations than variables.
- Fixed error from `solve.QP` when there are large values in the response.
- Now `checkCmat` also checks if there are "zero" constraints.
- Fixed a bug in `simulCoef`. Now includes a switch to simulate under the constrained or unconstrained model.

# 0.3.0

## New features
- Added `check_cmat` and `coef_simu` to the list of exported functions as they can be useful for specific use cases.
- Added full documentation.

## Bug fixes
- In `check_cmat`, removed the call to `limSolve::nnls()` to be replaced by `coneproj::coneB` (also a NNLS solver) to reduce the number of dependencies.

# 0.2.1

## New features
- Initialization of a short documentation for several functions.

## Bug fixes
- `cirls.control` now checks for constraint matrix irreducibility.
- `cirls.control` is now exported.
- changed the function to determine redundant constraints in `check_cmat`
- No warning message on row rank of `Cmat` anymore
- `Cmat` is now checked only when there are more than one row

# 0.2.0

## New features
- Method vcov for cirls to compute corrected covariance matrices
- Method confint for cirls to compute feasible confidence intervals
- Added several QP solvers: quadprog (the original one), osqp and coneproj.

## Changes
- A warning is now displayed when Cmat is not of full row rank
- vcov and confint return NA matrices if Cmat is not of full row rank
- Changed residual df computation to account for active constraints
- Replaced bvec by lb (lower bound) and ub (upper bound). Allows equality constraints.
- Added cirls class to glm output

## Bug fixes
- cirls.fit has the same behaviour as glm.fit when model is singular: fill with NAs
