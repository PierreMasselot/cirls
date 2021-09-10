
# Cmat is either a matrix with the right number of columns, or a named list of 
#   matrices for specific terms
cirls.control <- function (epsilon = 1e-08, maxit = 25, trace = FALSE,
  Cmat = NULL, bvec = NULL) 
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
    if (is.list(Cmat)){
      Cmat <- clist2cmat(Cmat)
    } else {
      if (ncol(Cmat) != ncol(get("x", envir = parent.frame()))){
        stop("Cmat must have the same number of columns as the design matrix")
      }
    }
  }
  if (is.null(bvec))
    bvec <- rep(0.0001, nrow(Cmat))
  if (NROW(bvec) != nrow(Cmat)){
    bvec <- rep_len(bvec, nrow(Cmat))
    warning("bvec has been recycled because of incompatibility with Cmat")
  }
  zerocons <- rowSums(Cmat != 0) == 0
  if (any(zerocons)){
    warning("constraints containing only zeros have been removed")
    Cmat <- Cmat[!zerocons,]
    bvec <- bvec[!zerocons]
  }
  list(epsilon = epsilon, maxit = maxit, trace = trace, Cmat = Cmat, 
    bvec = bvec)
}
