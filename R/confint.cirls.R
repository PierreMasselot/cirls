################################################################################
#
#   Confidence intervals
#
################################################################################

#' Simulate coefficients, calculate Confidence Intervals and Variance-Covariance Matrix for a `cirls` object.
#'
#' @description
#' `confint` computes confidence intervals for one of more parameters in a GLM fitted via [cirls.fit][cirls.fit()]. `vcov` compute the variance-covariance matrix of the parameters. Both methods are based on `coef_simu` that simulates coefficients from a Truncated Multivariate Normal distribution. These methods supersede the default [confint][stats::confint()] and [vcov][stats::vcov()] methods for `cirls` objects.
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
#' Both methods rely on the fact that \eqn{C\hat{\beta}} (with \eqn{C} the constraint matrix) follows a *Truncated Multivariate Normal* distribution
#' \deqn{C\hat{\beta} \sim TMVN(C\beta, CVC^T), l, u}
#' where TMVN represents a truncated Multivariate Normal distribution. \eqn{C} is the constraint matrix (`object$control$Cmat`) with bound \eqn{l} and \eqn{u}, while \eqn{V} is the unconstrained variance-covariance matrix (such as returned by `vcov.glm`).
#'
#' `coef_simu` simulates from the TMVN above and transforms back the realisations into the coefficients space. These realisations are then used by the `confint` and `vcov` methods which compute empirical quantiles and variance-covariance matrix, respectively. `coef_simu` is called internally by `confint` and `vcov` and doesn't need to be used directly, but it can be used to check other summaries of the coefficients distribution.
#'
#' @note
#' These methods only work when `Cmat` is of full row rank. If not the case, `Cmat` can be inspected through [check_cmat()].
#'
#' @returns
#' For `confint`, a two-column matrix with columns giving lower and upper confidence limits for each parameter.
#'
#' For `vcov`, a matrix of the estimated covariances between the parameter estimates of the model.
#'
#' For `coef_simu`, a matrix with `nsim` rows containing simulated coefficients.
#'
#' @references
#' Geweke, J.F., 1996. Bayesian Inference for Linear Models Subject to Linear Inequality Constraints, in: Lee, J.C., Johnson, W.O., Zellner, A. (Eds.), Modelling and Prediction Honoring Seymour Geisser. *Springer, New York, NY*, pp. 248–263. \doi{10.1007/978-1-4612-2414-3_15}
#'
#' Botev, Z.I., 2017, The normal law under linear restrictions: simulation and estimation via minimax tilting, *Journal of the Royal Statistical Society, Series B*, **79** (**1**), pp. 1–24.
#'
#' @seealso [rtmvnorm][TruncatedNormal::tmvnorm()] for the underlying routine to simulate from a TMVN. [check_cmat()] to check if the contraint matrix can be reduced.
#'
#' @example man/examples/new_methods.R
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
