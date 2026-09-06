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

  # Remove if aliased
  aliased <- stats::summary.glm(object)$aliased
  if (!complete) parm <- parm[!aliased[parm]]

  # simulate from truncated multivariate normal
  seed <- if ("seed" %in% names(dots)) dots$seed else NULL
  simures <- simulCoef(object, nsim = nsim, complete = TRUE, seed = seed)

  # Compute limits
  confint.sim.cirls(simures, parm = parm, level = level)
}
