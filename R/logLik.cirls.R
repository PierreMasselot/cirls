################################################################################
#
#  log-Likelihoos
#
################################################################################

#' Extract Log-Likelihood
#'
#' @description
#' Extracts the log-likelihood of a fitted `cirls` object with associated degrees of freedom. Typically used with [AIC][stats::AIC()].
#'
#' @param object A `cirls` object.
#' @param df Character. The type of degrees of freedom to assign to the log-Likelihood. Default to expected degrees of freedom `"edf"`. See [edf][edf()] for a description of possible degrees of freedom.
#' @param ... Additional arguments to be passed to [edf][edf()] to compute degrees of freedom.
#'
#' @details
#' The argument `df` provide the type of degrees of freedom attributed to the returned log-likelihood value. This is typically used in the computation of [AIC][stats::AIC()] and [BIC][stats::BIC()] and changing the degrees of freedom can ultimately change the values of the information criteria. By default, the expected number of freedom given the constraints is used. See [edf][edf()] for details on the computation and for the returned types of degrees of freedom.
#'
#' @returns
#' A numeric value of class `logLik` with attributes `df` (degrees of freedom, see details) and `nobs` (number of observations used in the estimation).
#'
#' @seealso [edf][edf()] to compute expected degrees of freedom.
#'
#' @example inst/examples/ex_warming_logLik.R
#'
#' @export
logLik.cirls <- function(object, df = "edf", ...){

  # Extract dfs
  dfvec <- edf(object, ...)

  # Check the type of df works
  df <- match.arg(df, names(dfvec))

  # Compute logLik
  p <- dfvec["odf"]
  val <- p - object$aic/2

  # Compute expected reduction in df due to constraints
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- dfvec[df]
  class(val) <- "logLik"
  val
}
