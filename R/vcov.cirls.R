################################################################################
#
#    Variance-covariance matrix
#
################################################################################

#' @rdname confint.cirls
#' @export
vcov.cirls <- function(object, complete = TRUE, nsim = 1000, ...)
{
  aliased <- summary(object)$aliased

  # Check constraint matrix
  Cmat <- object$Cmat
  rowrk <- qr(t(Cmat))$rank
  if (nrow(Cmat) > rowrk){
    warning("Cannot perform inference because Cmat is not full row rank")
    v <- matrix(NA, length(aliased), length(aliased),
      dimnames = list(names(aliased), names(aliased)))
    if (!complete) v <- v[which(!aliased), which(!aliased)]
    return(v)
  }

  # Simulate from truncated multivariate normal
  simures <- coef_simu(object, nsim = nsim, complete = complete)

  # Compute empirical variance
  v <- stats::var(simures)

  # Return
  return(v)
}

