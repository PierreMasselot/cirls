################################################################################
#
# Shape constraint matrix method:
# ns method
#
################################################################################

#' @export
cons_shape.ns <- function(x, diff = 0, sign = 1, ...){

  # Create the B-spline constraint matrix
  Cmat <- cons_shape.bs(x, diff, sign)

  # Adjust the boundary bases
  ord <- attr(x, "degree") + 1
  knots <- c(rep(attr(x, "Boundary.knots"), ord), attr(x, "knots"))
  knots <- sort(knots)
  const <- splines::splineDesign(knots, attr(x, "Boundary.knots"), ord = ord,
    derivs = c(ord, ord) / 2)
  if (!attr(x, "intercept")) const <- const[,-1]
  qr.const <- qr(t(const))
  Cmat <- as.matrix((t(qr.qty(qr.const, t(Cmat))))[, -(1L:2L)])

  # Return
  chkc <- check_cmat(Cmat)
  if (length(chkc$redundant) > 0) Cmat <- Cmat[-chkc$redundant, , drop = F]
  Cmat
}
