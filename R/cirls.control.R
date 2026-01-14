#' Parameters controlling CIRLS fitting
#'
#' @description Function controlling the [cirls.fit][cirls.fit()] algorithm. Typically only used internally with arguments passed to `...` in [glm][glm()], but may be used to construct a `control` argument to either function.
#'
#' @param constr A formula specifying constraints to be applied to specific terms in the model. See details in [buildCmat][buildCmat()] for constraint specification.
#' @param Cmat Constraint matrix specifying the linear constraints applied to coefficients. Can also be provided as a list of matrices for specific terms. See details in [buildCmat][buildCmat()] for constraint specification.
#' @param lb,ub Lower and upper bound vectors for the linear constraints. If not provided, defaults to 0 and Inf, respectively.
#' @param epsilon Positive convergence tolerance. The algorithm converges when the relative change in deviance is smaller than `epsilon`.
#' @param maxit Integer giving the maximal number of CIRLS iterations.
#' @param trace Logical indicating if output should be produced for each iteration.
#' @param qp_solver The quadratic programming solver. One of `"quadprog"` (the default), `"osqp"` or `"coneproj"`.
#' @param qp_pars List of parameters specific to the quadratic programming solver. See respective packages help.
#'
#' @details
#' The `control` argument of [glm][stats::glm()] is by default passed to the `control` argument of [cirls.fit][cirls.fit()], which uses its elements as arguments for [cirls.control][cirls.control()]: the latter provides defaults and sanity checking. The control parameters can alternatively be passed through the `...` argument of [glm][stats::glm()].
#'
#' ## Constraint specification
#'
#' Constraint specification through the `constr`, `Cmat`, `lb` and `ub` argument is fully detailed in the help of the [buildCmat][buildCmat()] functions.
#'
#' ## Quadratic programming solvers
#'
#' The function [cirls.fit][cirls.fit()] relies on a quadratic programming solver. Several solver are currently available.
#' - `"quadprog"` (the default) performs a dual algorithm to solve the quadratic program. It relies on the function [solve.QP][quadprog::solve.QP()].
#' - `"osqp"` solves the quadratic program via the Alternating Direction Method of Multipliers (ADMM). Internally calls the function [solve_osqp][osqp::solve_osqp()].
#' - `"coneproj"` solves the quadratic program by a cone projection method. It relies on the function [qprog][coneproj::qprog()].
#'
#' Each solver has specific parameters that can be controlled through the argument `qp_pars`. Sensible defaults are set within [cirls.control][cirls.control()] and the user typically doesn't need to provide custom parameters. `"quadprog"` is set as the default being generally more reliable than the other solvers. `"osqp"` is faster but can be less accurate, in which case it is recommended to increase convergence tolerance at the cost of speed.
#'
#' @returns A named list containing arguments to be used in [cirls.fit][cirls.fit()].
#'
#' @seealso Specification of a [cirls][cirls-package()] model and constraint specification in [buildCmat][buildCmat()].
#'
#' @example inst/examples/ex_london_nonneg.R
#'
#' @export
cirls.control <- function (constr = NULL, Cmat = NULL, lb = NULL, ub = NULL,
  epsilon = 1e-08, maxit = 25, trace = FALSE, qp_solver = "quadprog",
  qp_pars = list())
{

  # Check valid convergence parameters
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")

  #----- Check constraints

  # Check that Cmat/lb/ub are consistent
  cm <- list(Cmat, lb, ub)
  if (any(sapply(cm, is.numeric)) && any(sapply(cm, is.list)))
    stop("Cmat/lb/ub must be all either matrix/vector or named lists")

  # Check that constr is not used when Cmat/lb/ub passed as full
  if (any(sapply(cm, is.numeric)) && !is.null(constr))
    warning("constr will be ignored when Cmat/lb/ub are matrix/vector for full model")

  # Check `constr` is a formula or can be coerced as one
  if (!(inherits(constr, "formula") || is.null(constr))){
    stop("constr must be provided as a formula")
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
