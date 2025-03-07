#' Parameters controlling CIRLS fitting
#'
#' @description Internal function controlling the [glm][stats::glm()] fit with linear constraints. Typically only used internally by [cirls.fit][cirls.fit()], but may be used to construct a control argument.
#'
#' @param epsilon Positive convergence tolerance. The algorithm converges when the relative change in deviance is smaller than `epsilon`.
#' @param maxit Integer giving the maximal number of CIRLS iterations.
#' @param trace Logical indicating if output should be produced for each iteration.
#' @param Cmat Constraint matrix specifying the linear constraints applied to coefficients. Can also be provided as a list of matrices for specific terms.
#' @param lb,ub Lower and upper bound vectors for the linear constraints. Identical values in `lb` and `ub` identify equality constraints. As for `Cmat` can be provided as a list of terms. If some terms are provided in `Cmat` but not in `lb` or `ub`, default values of 0 and Inf will be used, respectively.
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
cirls.control <- function (epsilon = 1e-08, maxit = 25, trace = FALSE,
  Cmat = NULL, lb = 0L, ub = Inf, qp_solver = "quadprog", qp_pars = list())
{
  # Check valid convergence parameters
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")

  #----- Constraints preparation

  # Chech Cmat has been provided
  if (is.null(Cmat)){
    stop(paste0("Cmat must be provided. Note that 'control' and '...'",
      " cannot be used at the same time, with 'control' being prioritised."))
  }

  # Get objects that cannot be passed through arguments
  mt <- get("mt", envir = parent.frame(2))
  x <- get("x", envir = parent.frame())

  # If Cmat is a list, build matrix and bound vectors
  if (is.list(Cmat)){

    # Transform
    Cmat <- clist2cmat(Cmat, mt, x)
    ct <- attr(Cmat, "terms")

    # Match bounds
    bounds <- Map(function(b, def){
      if (is.list(b)){
        if (any(!names(b) %in% names(ct)))
          warning("'ub' and 'lb' terms not found in 'Cmat' are dropped")
        bvec <- rep(def, NROW(Cmat))
        for (i in seq_along(b)) bvec[ct[[names(b)[i]]]] <- b[[i]]
      } else {
        bvec <- rep_len(b, nrow(Cmat))
      }
      bvec
    }, b = list(lb, ub), def = c(0, Inf))
    lb <- bounds[[1]]
    ub <- bounds[[2]]
  } else {
    if (ncol(Cmat) != ncol(x)){
      stop("Cmat must have the same number of columns as the design matrix")
    }
    if (is.list(ub) || is.list(lb)) stop(paste0("when 'Cmat' ",
      "is provided as a matrix, 'lb' and 'ub' cannot be lists"))
    lb <- rep_len(lb, nrow(Cmat))
    ub <- rep_len(ub, nrow(Cmat))
  }

  # Check that bounds are well specified
  if (any(lb > ub)){
    warning("lb should be lower than (or equal) ub")
    lb2 <- pmin(lb, ub)
    ub2 <- pmax(lb, ub)
    lb <- lb2
    ub <- ub2
  }

  # Check irreducibility of Cmat
  if (nrow(Cmat) > 1){
    chkc <- checkCmat(Cmat)
    if (length(chkc$redundant) > 0){
      Cmat <- Cmat[-chkc$redundant,,drop = F]
      lb <- lb[-chkc$redundant]
      ub <- ub[-chkc$redundant]
      warning(paste0("Redundant constraints removed from Cmat: ",
        paste(chkc$redundant, collapse = ", ")))
    }
    if (length(chkc$equality) > 0){
      warning(paste0("Underlying equality constraints: ",
        paste(chkc$equality, collapse = ", "), ". ",
        "Consider using lb and ub to set equality constraints instead."))
    }
  }

  #----- Other parameters

  # Prepare QP solver
  qp_solver <- match.arg(qp_solver, c("quadprog", "osqp", "coneproj"))
  qp_pars <- do.call(sprintf("%s.def", qp_solver), qp_pars)

  # Return
  list(epsilon = epsilon, maxit = maxit, trace = trace, Cmat = Cmat,
    lb = lb, ub = ub, qp_solver = qp_solver, qp_pars = qp_pars)
}
