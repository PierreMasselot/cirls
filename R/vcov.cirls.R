################################################################################
#
#    Variance-covariance matrix
#
################################################################################

#' @rdname simulCoef
#' @order 3
#' @export
vcov.cirls <- function(object, complete = TRUE, nsim = 1000, trunc = TRUE,
  ...)
{
  aliased <- summary(object)$aliased
  dots <- list(...)

  if (trunc){

    # simulate from truncated multivariate normal
    seed <- if ("seed" %in% names(dots)) dots$seed else NULL
    simures <- simulCoef(object, nsim = nsim, complete = complete,
      seed = seed, constrained = TRUE)

    # Compute empirical variance
    v <- stats::var(simures)
  } else {

    # If trunc = FALSE return the unmodified vcov
    v <- stats::.vcov.aliased(aliased, stats::summary.glm(object)$cov.scaled,
      complete = complete)

  }
  # Return
  return(v)
}

