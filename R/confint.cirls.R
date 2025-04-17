################################################################################
#
#   Confidence intervals
#
################################################################################

#' @rdname simulCoef
#' @order 2
#' @export
confint.cirls <- function(object, parm, level = 0.95, nsim = 1000,
  complete = TRUE, ...)
{

  dots <- list(...)

  # Select coefficients
  aliased <- stats::summary.glm(object)$aliased
  pnames <- names(aliased)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]

  # Remove if aliased
  if (!complete) parm <- parm[!aliased[parm]]

  # simulate from truncated multivariate normal
  seed <- if ("seed" %in% names(dots)) dots$seed else NULL
  simures <- simulCoef(object, nsim = nsim, complete = TRUE, seed = seed)

  # Compute limits
  lims <- c((1 - level) / 2, level + (1 - level) / 2)
  res <- t(apply(simures[, parm, drop = F], 2, stats::quantile, lims,
    na.rm = TRUE))
  colnames(res) <- c("low", "high")

  # Return
  return(res)
}
