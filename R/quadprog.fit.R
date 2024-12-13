################################################################################
#
# Function to fit model by quadprog
#
################################################################################

quadprog.fit <- function(Dmat, dvec, Cmat, lb, ub, qp_pars){
  #----- Construct Cmat and bvec from lb and ub
  # Get equality constraints
  iseq <- lb == ub
  meq <- sum(iseq)
  Amat <- Cmat[iseq,]
  bvec <- lb[iseq]
  cmap <- which(iseq)

  # Get lb constraints
  lbcons <- lb[!iseq] > -Inf
  Amat <- rbind(Amat, Cmat[!iseq,,drop = F][lbcons,])
  bvec <- c(bvec, lb[!iseq][lbcons])
  cmap <- c(cmap, which(!iseq)[lbcons])

  # Get ub constraints
  ubcons <- ub[!iseq] < Inf
  Amat <- rbind(Amat, -Cmat[!iseq,,drop = F][ubcons,])
  bvec <- c(bvec, -ub[!iseq][ubcons])
  cmap <- c(cmap, which(!iseq)[ubcons])

  #----- Fit QP

  # Normalise the matrices to avoid erors due to huge numbers
  sc <- norm(Dmat, "2")
  Dmat <- Dmat / sc
  dvec <- dvec / sc

  # Fit
  res <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq)

  # Extract active constraints
  iact <- cmap[res$iact]

  # Return
  list(solution = res$solution, iterations = res$iterations[1],
    iact = unique(iact))
}
