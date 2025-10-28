################################################################################
#
# Shape constraint matrix method:
# ns method
#
################################################################################

#' @export
shapeConstr.ns <- function(x, shape, ...){

  # Create the B-spline constraint matrix
  cm <- shapeConstr.bs(x, shape, ...)

  # Adjust the boundary bases
  ord <- attr(x, "degree") + 1
  knots <- c(rep(attr(x, "Boundary.knots"), ord), attr(x, "knots"))
  knots <- sort(knots)
  const <- splines::splineDesign(knots, attr(x, "Boundary.knots"), ord = ord,
    derivs = c(2, 2))
  if (!attr(x, "intercept")) const <- const[,-1]
  qr.const <- qr(t(const))
  Cmat <- as.matrix((t(qr.qty(qr.const, t(cm$Cmat))))[, -(1L:2L)])

  # Constraining of NS can create some redundant constraints
  Cmat <- Cmat[!checkCmat(Cmat)$redundant,, drop = FALSE]

  # Bounds
  lb <- rep(0, NROW(Cmat))
  ub <- rep(Inf, NROW(Cmat))

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
