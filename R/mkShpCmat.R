################################################################################
#
# Make Constraint matrix for shape-constrained splines
#
################################################################################

mkShpCmat <- function(shape, knots, ord, intercept,
  mid, rng = c(-Inf, Inf), msg = FALSE)
{

  #----- Check parameters

  # Check shape and extract sign and diff parameters
  cpars <- chkshp(shape, ord)

  # Check range parameter
  rng <- chkrng(rng, knots, msg = msg)

  #----- Create constraint matrix

  # Create submatrix for each shape
  Cmat <- lapply(cpars, function(cp) mkDmat(cp[1], cp[2], knots = knots,
    ord = ord, lower = rng[1], upper = rng[2]))

  # Put together
  Cmat <- do.call(rbind, Cmat)

  # Optionally remove intercept
  if (!intercept) Cmat <- Cmat[,-1, drop = F]

  # Remove redundant constraints
  Cmat <- reduceCons(Cmat, warn = FALSE)$Cmat

  #----- Return

  # Add bounds attributes
  clist <- list(Cmat = Cmat, lb = rep(0, NROW(Cmat)), ub = rep(Inf, NROW(Cmat)))

  # Return
  clist
}
