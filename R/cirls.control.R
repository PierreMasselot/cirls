
# Cmat is either a matrix with the right number of columns, or a named list of
#   matrices for specific terms
cirls.control <- function (epsilon = 1e-08, maxit = 25, trace = FALSE,
  Cmat = NULL, lb = 0L, ub = Inf, qp_solver = "osqp", qp_pars = list())
{
  # Check valid convergence parameters
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  # Chech Cmat and prepare it
  if (is.null(Cmat)){
    stop("Cmat must be provided")
  } else {
    # Get objects that cannot be passed through arguments
    mt <- get("mt", envir = parent.frame(2))
    x <- get("x", envir = parent.frame())
    if (is.list(Cmat)){
      Cmat <- clist2cmat(Cmat, mt, x)
    } else {
      if (ncol(Cmat) != ncol(x)){
        stop("Cmat must have the same number of columns as the design matrix")
      }
    }
  }
  # Check lb and recycle if needed
  if (NROW(lb) != nrow(Cmat)){
    lb <- rep_len(lb, nrow(Cmat))
  }
  # Same with ub and recycle if needed
  if (NROW(ub) != nrow(Cmat)){
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
  # Check that all constraints are proper
  zerocons <- rowSums(Cmat != 0) == 0
  if (any(zerocons)){
    warning("constraints containing only zeros have been removed")
    Cmat <- Cmat[!zerocons,]
    lb <- lb[!zerocons]
    ub <- ub[!zerocons]
  }
  # Check row rank of Cmat
  rowrk <- qr(t(Cmat))$rank
  if (nrow(Cmat) > rowrk){
    warning(paste0("Cmat does not have full row rank and inference won't ",
      "be possible. Check for possibly redundant constraints"))
  }
  # Prepapre QP solver
  qp_solver <- match.arg(qp_solver, c("quadprog", "osqp", "coneproj"))
  qp_pars <- do.call(sprintf("%s.def", qp_solver), qp_pars)
  # Return
  list(epsilon = epsilon, maxit = maxit, trace = trace, Cmat = Cmat,
    lb = lb, ub = ub, qp_solver = qp_solver, qp_pars = qp_pars)
}
