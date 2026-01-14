################################################################################
#
# Shape constraint matrix method:
# Default method
#
################################################################################

#' @rdname shapeConstr
#' @export
shapeConstr.factor <- function(x, shape, intercept = FALSE, ...) {

  # Levels
  ord <- nlevels(x)

  # Check parameters
  cpars <- chkshp(shape, ord)

  # Create constraint matrices
  knots <- seq_len(2 * ord)
  Cmat <- lapply(cpars, function(cp) dmat(cp[1], cp[2], knots, ord))
  Cmat <- do.call(rbind, Cmat)

  # Check if the "intercept" is not included, contrasts are applied
  # NB: no contrast is applied in R when the model does not include an intercept
  if (isFALSE(intercept)){
    ctr <- stats::contrasts(x)
    Cmat <- Cmat %*% ctr
  }

  # Remove redundant constraints
  Cmat <- Cmat[!checkCmat(Cmat, warn = FALSE)$redundant,, drop = FALSE]

  # Bounds
  lb <- rep(0, NROW(Cmat))
  ub <- rep(Inf, NROW(Cmat))

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
