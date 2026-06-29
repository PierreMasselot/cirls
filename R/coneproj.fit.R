################################################################################
#
# Function to fit model by coneproj
#
################################################################################

coneproj.fit <- function(Rmat, effects, Cmat, lb, ub, qp_pars){

  #----- Construct Cmat and bvec from lb and ub
  if (NROW(Cmat) > 0){
    # Get lb constraints
    lbcons <- lb > -Inf
    Amat <- Cmat[lbcons,,drop = F]
    bvec <- lb[lbcons]
    cmap <- which(lbcons)

    # Get ub constraints
    ubcons <- ub < Inf
    Amat <- rbind(Amat, -Cmat[ubcons,])
    bvec <- c(bvec, -ub[ubcons])
    cmap <- c(cmap, which(ubcons))
  } else {
    stop("coneproj does not accept empty Cmat, use another solver")
  }

  #----- Fit

  # Compute matrices
  qmat <- crossprod(Rmat)
  cvec <- crossprod(effects, Rmat)

  # Fit
  res <- coneproj::qprog(q = qmat, c = cvec, amat = Amat, b = bvec,
    msg = qp_pars$msg)

  # Get active constraints
  # iact <- which(mapply(function(x, y) isTRUE(all.equal(x, y)),
  #   Amat %*% res$thetahat - bvec, 0))
  iact <- res$face

  # Return
  list(solution = res$thetahat, iterations = res$steps,
    iact = unique(cmap[iact]))
}
