################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# ns method
#
################################################################################

# By default deg is set to the degree of the spline, for smoother convergence

#' @rdname boundConstr.bs
#' @order 2
#' @export
boundConstr.ns <- function(x, ...)
{

  attrs <- attributes(x)[c("degree", "knots", "Boundary.knots", "intercept")]
  attrs$x <- seq(attrs$Boundary.knots[1], attrs$Boundary.knots[2],
    length.out = 10)
  xbs <- do.call(splines::bs, attrs)

  # Constraint matrix
  cm <- boundConstr(xbs, ...)

  # Adjust the constraint matrix with ns bound conditions
  ord <- attr(x, "degree") + 1
  knots <- c(rep(attr(x, "Boundary.knots"), ord), attr(x, "knots"))
  knots <- sort(knots)
  const <- splines::splineDesign(knots, attr(x, "Boundary.knots"), ord = ord,
    derivs = c(2, 2))
  if (!attr(x, "intercept")) const <- const[,-1]
  qr.const <- qr(t(const))
  Cmat <- as.matrix((t(qr.qty(qr.const, t(cm$Cmat))))[, -(1L:2L), drop = F])

  # Constraining of NS can create some redundant constraints
  chk <- reduceCons(Cmat = Cmat, lb = cm$lb, ub = cm$ub, warn = FALSE)

  # Return
  chk[c("Cmat", "lb", "ub")]
}
