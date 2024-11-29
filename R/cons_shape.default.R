################################################################################
#
# Shape constraint matrix method:
# Default method
#
################################################################################

#' @export
cons_shape.default <- function(x, diff = 0, sign = 1, intercept = F, ...) {

  # Number of bases
  ord <- ncol(x) + !intercept

  # Check parameters
  cpars <- chkcpars(diff, sign, intercept, ord)
  ncons <- length(cpars$diff)

  # Initialise diagonal matrices
  Cmat <- replicate(ncons, diag(ord), simplify = FALSE)

  # Compute the difference when diff > 0
  ind <- cpars$diff > 0
  Cmat[ind] <- Map(function(C, s, d) s * diff(C, diff = d),
    Cmat[ind], cpars$sign[ind], cpars$diff[ind])

  # Put together and remove redundant constraints
  Cmat <- do.call(rbind, Cmat)
  if (!intercept) Cmat <- Cmat[,-1]
  chkc <- check_cmat(Cmat)
  if (length(chkc$redundant) > 0) Cmat <- Cmat[-chkc$redundant, , drop = F]

  # Return
  Cmat
}
