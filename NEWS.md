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

## Bug correction
- cirls.fit has the same behaviour as glm.fit when model overdetermined: fill with NAs
