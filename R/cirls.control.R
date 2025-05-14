#' Parameters controlling CIRLS fitting
#'
#' @description Internal function controlling the [glm][stats::glm()] fit with linear constraints. Typically only used internally by [cirls.fit][cirls.fit()], but may be used to construct a control argument.
#'
#' @param constr A formula specifying constraints to be applied to specific terms in the model.
#' @param Cmat Constraint matrix specifying the linear constraints applied to coefficients. Can also be provided as a list of matrices for specific terms.
#' @param lb,ub Lower and upper bound vectors for the linear constraints. Identical values in `lb` and `ub` identify equality constraints. As for `Cmat` can be provided as a list of terms. If some terms are provided in `Cmat` but not in `lb` or `ub`, default values of 0 and Inf will be used, respectively.
#' @param epsilon Positive convergence tolerance. The algorithm converges when the relative change in deviance is smaller than `epsilon`.
#' @param maxit Integer giving the maximal number of CIRLS iterations.
#' @param trace Logical indicating if output should be produced for each iteration.
#' @param qp_solver The quadratic programming solver. One of `"quadprog"` (the default), `"osqp"` or `"coneproj"`.
#' @param qp_pars List of parameters specific to the quadratic programming solver. See respective packages help.
#'
#' @details
#' The `control` argument of [glm][stats::glm()] is by default passed to the `control` argument of [cirls.fit][cirls.fit()], which uses its elements as arguments for [cirls.control][cirls.control()]: the latter provides defaults and sanity checking. The control parameters can alternatively be passed through the `...` argument of [glm][stats::glm()]. See [glm.control][stats::glm.control()] for details on general GLM fitting control, and [cirls.fit][cirls.fit()] for details on arguments specific to constrained GLMs.
#'
#' @returns A named list containing arguments to be used in [cirls.fit][cirls.fit()].
#'
#' @seealso the main function [cirls.fit][cirls.fit()], and [glm.control][stats::glm.control()].
#'
#' @example man/examples/cirls.control.R
#'
#' @export
cirls.control <- function (constr = NULL, Cmat = NULL, lb = 0L, ub = Inf,
  epsilon = 1e-08, maxit = 25, trace = FALSE,
  qp_solver = "quadprog", qp_pars = list())
{
  # Check valid convergence parameters
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")

  #----- Check constraints

  # Check that bounds are well specified
  if (any(lb > ub)){
    warning("lb cannot greater than ub")
    lb2 <- pmin(lb, ub)
    ub2 <- pmax(lb, ub)
    lb <- lb2
    ub <- ub2
  }

  # Check `constr` is a formula or can be coerced as one
  if (!(inherits(constr, "formula") || is.null(constr))){
    stop("constr must be provided as a formula")
  }

  # Check Cmat is either a matrix or a list of numerical matrix/vectors
  if (!(is.matrix(Cmat) || is.null(Cmat))){
    if (!all(sapply(Cmat, is.numeric))) stop(
      "Cmat must a numeric matrix or a list of numeric matrices")
  }

  #----- Other parameters

  # Prepare QP solver
  qp_solver <- match.arg(qp_solver, c("quadprog", "osqp", "coneproj"))
  qp_pars <- do.call(sprintf("%s.def", qp_solver), qp_pars)

  # Return
  list(constr = constr, Cmat = Cmat, lb = lb, ub = ub,
    epsilon = epsilon, maxit = maxit, trace = trace,
    qp_solver = qp_solver, qp_pars = qp_pars)
}
