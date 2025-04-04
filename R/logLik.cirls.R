################################################################################
#
#  log-Likelihoos
#
################################################################################

#' Log-Likelihood for a fitted `cirls` object
#'
#' @description
#' Extracts the log-likelihood for a fitted `cirls` object to be typically used by [AIC][stats::AIC()].
#'
#' @param object A `cirls` object.
#' @param df The type of degrees of freedom to assign to the log-Likelihood. Default to expected degrees of freedom. See [edf()].
#' @param ... Arguments to be passed to [edf][edf()] to compute degrees of freedom.
#'
#' @details
#' The argument `df` provide the type of degrees of freedom attributed to the returned log-likelihood value. This is typically used in the computation of [AIC][stats::AIC()] and [BIC][stats::BIC()] and changing the degrees of freedom can ultimately change the values of the information criteria. By default, the expected number of freedom given the constraints is used. See [edf][edf()] for details on the computation and for the returned types of degrees of freedom.
#'
#' @returns
#' A numeric value of class `logLik` with attributes `df` (degrees of freedom, see details) and `nobs` (number of observations used in the estimation).
#'
#' @seealso [edf][edf()] to compute expected degrees of freedom.
#'
#' @export
logLik.cirls <- function(object, df = "edf", ...){

  df <- match.arg(df, c("edf", "odf", "udf"))

  # Extract dfs
  dfvec <- edf(object, ...)

  # Compute logLik
  p <- dfvec["edf"]
  val <- p - object$aic/2

  # Compute expected reduction in df due to constraints
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- dfvec[df]
  class(val) <- "logLik"
  val
}
