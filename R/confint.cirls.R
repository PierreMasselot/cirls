################################################################################
#
#   Confidence intervals
#
################################################################################

#' Calculate Confidence Intervals and Variance-Covariance Matrix for a `cirls` object.
#'
#' @description
#' `confint` computes confidence intervals for one of more parameters in a GLM fitted via [cirls.fit][cirls.fit()]. `vcov` compute the variance-covariance matrix of the parameters. These methods supersede the default [confint][stats::confint()] and [vcov][stats::vcov()] methods for `cirls` objects.
#'
#' @param object A fitted `cirls` object.
#' @param parm A specification of which parameters to compute the confidence intervals for. Either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required.
#' @param nsim The number of simulations to consider. Corresponds to `n` in [rtmvnorm][TruncatedNormal::tmvnorm()]. See details().
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @details
#' These functions are custom methods for [cirls][cirls.fit()] objects to supersede the default methods used for [glm][stats::glm()] objects.
#'
#' Both methods simulate `nsim` realisations from a Truncated Normal distribution in which the truncation is defined by the constraints. Variance-Covariance matrix and Confidence Intervals are then computed as the empirical variance-covariance matrix and quantiles, respectively.
#'
#' @returns
#' For `confint`, a two-column matrix with columns giving lower and upper confidence limits for each parameter.
#'
#' For `vcov`, a matrix of the estimated covariances between the parameter estimates of the model.
#'
#' @references
#' Geweke, J.F., 1996. Bayesian Inference for Linear Models Subject to Linear Inequality Constraints, in: Lee, J.C., Johnson, W.O., Zellner, A. (Eds.), Modelling and Prediction Honoring Seymour Geisser. *Springer, New York, NY*, pp. 248–263. [10.1007/978-1-4612-2414-3_15](https://doi.org/10.1007/978-1-4612-2414-3_15)
#'
#' Botev, Z.I., 2017, The normal law under linear restrictions: simulation and estimation via minimax tilting, *Journal of the Royal Statistical Society, Series B*, **79** (**1**), pp. 1–24.
#'
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
