################################################################################
#
#   Confidence intervals
#
################################################################################

#' @export
confint.cirls <- function(object, parm, level = 0.95, nsim = 1000, ...)
{

  # Select coefficients
  pnames <- names(stats::coef(object))
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]

  # Remove if aliased
  aliased <- stats::summary.glm(object)$aliased
  parma <- parm[!aliased]

  # Check constraint matrix
  Cmat <- object$Cmat
  rowrk <- qr(t(Cmat))$rank
  if (nrow(Cmat) > rowrk){
    warning(paste0("Cannot perform inference because Cmat is not full row rank. ",
      "Check for possibly redundant constraints"))
    return(matrix(NA, length(parm), 2, dimnames = list(parm, c("low", "high"))))
  }

  # Simulate from truncated multivariate normal
  simures <- coef_simu(object, nsim)

  # Compute limits
  lims <- c((1 - level) / 2, level + (1 - level) / 2)
  res <- t(apply(simures[, parma, drop = F], 2, stats::quantile, lims))
  colnames(res) <- c("low", "high")

  # Add aliased and return
  if(any(aliased)){
    ares <- matrix(NA, length(parm), 2, dimnames = list(parm, c("low", "high")))
    ares[parma,] <- res
    return(ares)
  } else {
    return(res)
  }
}
