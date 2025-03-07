################################################################################
#
#    Variance-covariance matrix
#
################################################################################

#' @rdname simulCoef
#' @order 3
#' @export
vcov.cirls <- function(object, complete = TRUE, nsim = 1000, constrained = TRUE,
  ...)
{
  aliased <- summary(object)$aliased
  dots <- list(...)

  if (constrained){

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

    # simulate from truncated multivariate normal
    seed <- if ("seed" %in% names(dots)) dots$seed else NULL
    simures <- simulCoef(object, nsim = nsim, complete = complete,
      seed = seed)

    # Compute empirical variance
    v <- stats::var(simures)
  } else {

    # If constrained = FALSE return the unmodified vcov
    v <- stats::.vcov.aliased(aliased, stats::summary.glm(object)$cov.scaled,
      complete = complete)

  }
  # Return
  return(v)
}

