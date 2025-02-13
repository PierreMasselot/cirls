################################################################################
#
# Shape constraint matrix method:
# Default method
#
################################################################################

#' @export
shapeConstr.default <- function(x, shape, intercept = FALSE, ...) {

  # Check parameters
  cpars <- chkshp(shape)

  # Create constraint matrices
  ord <- max(sapply(cpars, "[", 1)) + 1
  knots <- seq_len(ncol(x) + ord + !intercept)
  Cmat <- lapply(cpars, function(cp) dmat(cp[1], cp[2], knots, ord))
  Cmat <- do.call(rbind, Cmat)

  # Put together and remove redundant constraints
  if (!intercept) Cmat <- Cmat[,-1]
  chkc <- checkCmat(Cmat)
  if (length(chkc$redundant) > 0) Cmat <- Cmat[-chkc$redundant, , drop = F]

  # Return
  Cmat
}
