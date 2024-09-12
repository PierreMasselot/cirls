################################################################################
#
# Function to fit model by coneproj
#
################################################################################

coneproj.fit <- function(Dmat, dvec, Cmat, lb, ub, qp_pars){
  #----- Construct Cmat and bvec from lb and ub

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

  #----- Fit

  # Fit
  res <- coneproj::qprog(q = Dmat, c = dvec, amat = Amat, b = bvec,
    msg = qp_pars$msg)

  # Get active constraints
  # iact <- which(mapply(function(x, y) isTRUE(all.equal(x, y)),
  #   Amat %*% res$thetahat - bvec, 0))
  iact <- res$face

  # Return
  list(solution = res$thetahat, iterations = res$steps,
    iact = unique(cmap[iact]))
}
