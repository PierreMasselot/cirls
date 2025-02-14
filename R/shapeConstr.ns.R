################################################################################
#
# Shape constraint matrix method:
# ns method
#
################################################################################

#' @rdname shapeConstr
#' @export
shapeConstr.ns <- function(x, shape, ...){

  # Create the B-spline constraint matrix
  Cmat <- shapeConstr.bs(x, shape, ...)

  # Adjust the boundary bases
  ord <- attr(x, "degree") + 1
  knots <- c(rep(attr(x, "Boundary.knots"), ord), attr(x, "knots"))
  knots <- sort(knots)
  const <- splines::splineDesign(knots, attr(x, "Boundary.knots"), ord = ord,
    derivs = c(2, 2))
  if (!attr(x, "intercept")) const <- const[,-1]
  qr.const <- qr(t(const))
  Cmat <- as.matrix((t(qr.qty(qr.const, t(Cmat))))[, -(1L:2L)])

  # Return
  chkc <- checkCmat(Cmat)
  if (length(chkc$redundant) > 0) Cmat <- Cmat[-chkc$redundant, , drop = F]
  Cmat
}
