#' Unconstrained model
#'
#' @description
#' Takes a fitted `cirls` object and returns the corresponding unconstrained model. Works similarly to the [update][stats::update()] method by removing all argument specific to [cirls.fit][cirls.fit()] from the call and evaluate it.
#'
#' @param object Fitted 'cirls' object.
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

  #----- Extracting everything: cumbersome
  # # Extract the model.frame and control for glm
  # mf <- model.frame(object)
  # mt <- attr(mf, "terms")
  # control <- object$control[names(object$control) %in%
  #     methods::formalArgs(stats::glm.control)]
  #
  # # Fit with glm.fit
  # # NOTE: the starting values and singular.ok are ignored
  # fit <- glm.fit(x = model.matrix(mf), y = model.response(mf),
  #   weights = model.weights(mf), offset = model.offset(mf),
  #   family = object$family, control = control,
  #   intercept = attr(mt, "intercept") > 0L)
  #
  # # Put together as in glm
  # if (length(model.offset(mf)) && attr(mt, "intercept") > 0L) {
  #   fit2 <- glm.fit(x = model.matrix(mf)[, "(Intercept)", drop = FALSE],
  #     y = model.response(mf), mustart = fit$fitted.values,
  #     weights = model.weights(mf), offset = model.offset(mf),
  #     family = object$family, control = control, intercept = TRUE)
  #   if (!fit2$converged)
  #     warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
  #   fit$null.deviance <- fit2$deviance
  # }
  # if (!is.null(object$model)) fit$model <- mf
  # fit$na.action <- attr(mf, "na.action")
  # if (!is.null(object$x)) fit$x <- model.matrix(mf)
  # if (is.null(object$y)) fit$y <- NULL
  # structure(c(fit, list(call = object$call, formula = object$formula, terms = mt,
  #   data = object$data, offset = model.offset(mf), control = control,
  #   method = "glm.fit", contrasts = attr(model.matrix(mf), "contrasts"),
  #   xlevels = .getXlevels(mt, mf))), class = c(fit$class, c("glm", "lm")))

  #----- Modifying the call instead: works well with the data environment

  # Extract call and change the method argument
  call <- as.list(stats::getCall(object))

  # Remove arguments related to cirls only from the call
  cirlsargs <- c("method", setdiff(names(call)[-1],
    methods::formalArgs(stats::glm)))
  call[cirlsargs] <- NULL

  # Modify the control list if present
  if ("control" %in% names(call)){
    control <- as.list(call$control)
    remargs <- setdiff(names(control)[-1],
      methods::formalArgs(stats::glm.control))
    control[remargs] <- NULL
    call$control <- if (length(control) > 1) as.call(control) else NULL
  }

  # Evaluate
  eval(as.call(call), as.environment(object$data))
}
