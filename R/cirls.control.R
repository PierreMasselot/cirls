
# Cmat is either a matrix with the right number of columns, or a named list of
#   matrices for specific terms
cirls.control <- function (epsilon = 1e-08, maxit = 25, trace = FALSE,
  Cmat = NULL, bvec = 0L, bvectol = 1e-04)
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
  # Check bvec and recycle if needed
  if (NROW(bvec) != nrow(Cmat)){
    bvec <- rep_len(bvec, nrow(Cmat))
  }
  # tolerance level for bvec to avoid unfeasible solutions due to rounding
  #   errors
  bvec <- bvec + bvectol
  # Check that all constraints are proper
  zerocons <- rowSums(Cmat != 0) == 0
  if (any(zerocons)){
    warning("constraints containing only zeros have been removed")
    Cmat <- Cmat[!zerocons,]
    bvec <- bvec[!zerocons]
  }
  list(epsilon = epsilon, maxit = maxit, trace = trace, Cmat = Cmat,
    bvec = bvec, bvectol = bvectol)
}
