#' Unconstrained model
#'
#' @description
#' Takes a fitted `cirls` object and fits the corresponding unconstrained model.
#'
#' @param object A `cirls` object.
#'
#' @details
#' This function refits a `cirls` object, removing all constraints. This is primarly used by [simulCoef][simulCoef()] for inference, but it can also be used to easily compare a constrained and an unconstrained model.
#'
#' The function still fits the model through [cirls.fit][cirls.fit()], but using an empty constraint matrix. Therefore, it returns a `cirls` object that can use all the facilities provided by the `cirls` package. In this instance, the CIRLS algorithm reduces to a classical IRLS and the results are identical to a usual [glm][glm()] fitted with `glm.fit`.
#'
#' @note
#' If any starting values were provided to fit the `cirls` object, they are not transferred to the fitting of the unconstrained model.
#'
#' @returns A `cirls` object.
#'
#' @seealso [simulCoef][simulCoef()] to perform inference.
#'
#' @example inst/examples/ex_warming_factor.R
#'
#' @export
uncons <- function(object){

  #----- Extract the data and fit the model

  # Extract specific model components
  mf <- stats::model.frame(object)
  x <- stats::model.matrix(object)
  y <- object$y %||% stats::model.response(mf)
  control <- object$control
  control[["constr"]] <- NULL
  mt <- stats::terms(object) # Needed within cirls.fit
  intercept <- attr(mt, "intercept") > 0L

  # Fit with cirls.fit
  # ignores the warning message displayed when there is no constraint
  withCallingHandlers(
    fit <- cirls.fit(x = x, y = y,
      weights = object$prior.weights, etastart = object$etastart,
      offset = object$offset, family = object$family, control = control,
      intercept = intercept, singular.ok = object$singular.ok),
    warning = function(w) {
      if (startsWith(conditionMessage(w), "No constraint provided"))
        invokeRestart("muffleWarning")
    }
  )

  # Part to compute the null deviance in the case of an offset and no intercept
  if (length(object$offset) && attr(mt, "intercept") > 0L) {
    control2 <- utils::modifyList(control, list(Cmat = as.matrix(1)))
    fit2 <- cirls.fit(x = x[, "(Intercept)", drop = FALSE],
      y = y, mustart = fit$fitted.values,
      weights = object$prior.weights, offset = object$offset,
      family = object$family, control = control2, intercept = TRUE)
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }

  #----- Export everything

  # Extract call and change the lb and ub arguments
  call <- as.list(object$call)
  call[c("lb", "ub")] <- list(lb = -Inf, ub = Inf)

  # Also put the data as in the fitted cirls
  fit$model <- object[["model"]]
  fit$y <- object[["y"]]
  fit$x <- object[["x"]]

  # Final list
  structure(c(fit,
    list(call = as.call(call), formula = object$formula, terms = mt,
      data = object$data, offset = object$offset, control = object$control,
      method = "cirls.fit", contrasts = attr(x, "contrasts"),
      xlevels = stats::.getXlevels(mt, mf))),
    class = c(fit$class, c("glm", "lm")))
}
