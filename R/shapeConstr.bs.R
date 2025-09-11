################################################################################
#
# Shape constraint matrix method:
# bs method
#
################################################################################

#' @export
shapeConstr.bs <- function(x, shape, ...){

  # Extract specifications
  intercept <- attr(x, "intercept")
  ord <- attr(x, "degree") + 1
  knots <- c(rep(attr(x, "Boundary.knots"), ord), attr(x, "knots"))
  knots <- sort(knots)

  # Extract parameters
  cpars <- chkshp(shape, ord)

  # Create Constraint matrix
  Cmat <- lapply(cpars, function(cp) dmat(cp[1], cp[2], knots, ord))
  Cmat <- do.call(rbind, Cmat)

  # Remove redundant constraints
  if (!intercept) Cmat <- Cmat[,-1]
  Cmat <- Cmat[!checkCmat(Cmat)$redundant,, drop = FALSE]

  # Add bound attributes
  attr(Cmat, "lb") <- rep(0, NROW(Cmat))
  attr(Cmat, "ub") <- rep(Inf, NROW(Cmat))

  # Return
  Cmat
}
