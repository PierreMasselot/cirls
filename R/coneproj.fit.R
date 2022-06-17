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

  # get active constraints
  iact <- which((Amat %*% res$thetahat - bvec + 1) == 1)

  # Return
  list(solution = res$thetahat - 1 + 1, iterations = res$steps,
    iact = unique(cmap[iact]))
}
