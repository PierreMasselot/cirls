#' Unconstrained model
#'
#' @description
#' Takes a fitted `cirls` object and returns the corresponding unconstrained model. Works similarly to the [update][stats::update()] method by removing all argument specific to [cirls.fit][cirls.fit()] from the call and evaluate it.
#'
#' @param object Fitted 'cirls' object.
#'
#' @details
#'
#' ## Note on starting values
#'
#' If any starting values were provided to fit the `cirls` object, they are not transferred to the fitting of the unconstrained model.
#'
#' @returns A [glm][stats::glm()] object.
#'
#' @examples
#' # Generate some data
#' n <- 1000
#' betas <- c(0, 1, 2, -1, 1)
#' p <- length(betas)
#' x <- matrix(rnorm(n * p), n, p)
#' eta <- 5 + x %*% betas
#' y <- eta + rnorm(n, 0, .2)
#'
#' # Fit two cirls models, passing Cmat through the two different pathways
#' cinc <- diff(diag(p))
#' res1 <- glm(y ~ x, method = cirls.fit, Cmat = list(x = cinc))
#' res2 <- glm(y ~ x, method = cirls.fit, control = list(Cmat = list(x = cinc)))
#'
#' # 'Unconstrain' the models
#' uc1 <- uncons(res1)
#' uc2 <- uncons(res2)
#'
#' @export
uncons <- function(object){

  #----- Extract the data and fit the model

  # Extract the model.frame and control for glm
  mf <- model.frame(object)
  mt <- attr(mf, "terms")
  control <- object$control[names(object$control) %in%
      methods::formalArgs(stats::glm.control)]
  x <- model.matrix(mt, mf)
  y <- model.response(mf)
  weights <- model.weights(mf)
  offset <- model.offset(mf)

  # Fit with glm.fit
  # NOTE: the starting values and singular.ok are ignored
  fit <- glm.fit(x = x, y = y, weights = weights, offset = offset,
    family = object$family, control = control,
    intercept = attr(mt, "intercept") > 0L)

  # Part to compute the null deviance in the case of an offset and no intercept
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- glm.fit(x = x[, "(Intercept)", drop = FALSE],
      y = y, mustart = fit$fitted.values,
      weights = weights, offset = offset,
      family = object$family, control = control, intercept = TRUE)
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }

  #----- Modify the call to match the new model

  # Extract call and change the method argument
  call <- as.list(object$call)

  # Remove arguments related to cirls only from the call
  cirlsargs <- c("method", setdiff(names(call)[-1],
    methods::formalArgs(stats::glm)))
  call[cirlsargs] <- NULL

  # Modify the control list if present
  if ("control" %in% names(call)){
    clcontrol <- as.list(call$control)
    remargs <- setdiff(names(clcontrol)[-1],
      methods::formalArgs(stats::glm.control))
    clcontrol[remargs] <- NULL
    call$control <- if (length(clcontrol) > 1) as.call(clcontrol) else NULL
  }

  #----- Export everything

  # Also put the data as in the fitted cirls
  fit$model <- object[["model"]]
  fit$y <- object[["y"]]
  fit$x <- object[["x"]]

  # Final list
  structure(c(fit,
    list(call = as.call(call), formula = object$formula, terms = mt,
      data = object$data, offset = offset,
      control = do.call("glm.control", control),
      method = "glm.fit", contrasts = attr(x, "contrasts"),
      xlevels = .getXlevels(mt, mf))),
    class = c(fit$class, c("glm", "lm")))
}
