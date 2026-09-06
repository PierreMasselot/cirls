################################################################################
#
#   Confidence intervals for sim.cirls objects
#
################################################################################

#' @rdname simulCoef
#' @order 3
#' @export
confint.sim.cirls <- function(object, parm, level = 0.95, ...)
{

  # Select coefficients
  pnames <- colnames(object)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm)) {
    if (!attr(object, "complete"))
      warning(paste0("'object' has been simulated without",
        " aliased coefficients, so indices provided in 'parm' may not",
        " exactly correspond to the design matrix"))
    parm <- pnames[parm]
  }

  # Compute limits
  lims <- c((1 - level) / 2, level + (1 - level) / 2)
  res <- t(apply(object[, parm, drop = F], 2, stats::quantile, lims,
    na.rm = TRUE))
  colnames(res) <- c("low", "high")

  # Return
  res
}
