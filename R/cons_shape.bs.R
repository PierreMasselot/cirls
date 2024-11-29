################################################################################
#
# Shape constraint matrix method:
# bs method
#
################################################################################

#' @export
cons_shape.bs <- function(x, diff = 0, sign = 1, ...){

  # Extract specifications
  intercept <- attr(x, "intercept")
  ord <- attr(x, "degree") + 1
  knots <- c(rep(attr(x, "Boundary.knots"), ord), attr(x, "knots"))
  knots <- sort(knots)

  # Check parameters
  cpars <- chkcpars(diff, sign, intercept, ord)
  ncons <- length(cpars$diff)

  # Initialise diagonal matrices
  Cmat <- replicate(ncons, diag(length(knots) - ord), simplify = FALSE)

  # Create Constraint matrices
  ind <- cpars$diff > 0
  Cmat[ind] <- Map(function(s, d) s * dmat(d, knots, ord),
    cpars$sign[ind], cpars$diff[ind])

  # Put together and remove redundant constraints
  Cmat <- do.call(rbind, Cmat)
  if (!intercept) Cmat <- Cmat[,-1]
  chkc <- check_cmat(Cmat)
  if (length(chkc$redundant) > 0) Cmat <- Cmat[-chkc$redundant, , drop = F]

  # Return
  Cmat
}
