################################################################################
#
#    Variance-covariance matrix
#
################################################################################

#' @rdname confint.cirls
#' @export
vcov.cirls <- function(object, nsim = 1000, ...)
{
  aliased <- summary(object)$aliased

  # Check constraint matrix
  Cmat <- object$Cmat
  rowrk <- qr(t(Cmat))$rank
  if (nrow(Cmat) > rowrk){
    warning(paste0("Cannot perform inference because Cmat is not full row rank. ",
      "Check for possibly redundant constraints"))
    return(matrix(NA, length(aliased), length(aliased),
      dimnames = list(names(aliased), names(aliased))))
  }

  # Simulate from truncated multivariate normal
  simures <- coef_simu(object, nsim)

  # Compute empirical variance
  v <- stats::var(simures)

  # Add NAs for aliased coefficients
  if (any(aliased)){
    va <- matrix(NA, length(aliased), length(aliased),
      dimnames = list(names(aliased), names(aliased)))
    va[!aliased, !aliased] <- v
    return(va)
  } else{
    return(v)
  }
}

