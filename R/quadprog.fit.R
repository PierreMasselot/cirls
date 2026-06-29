################################################################################
#
# Function to fit model by quadprog
#
################################################################################

quadprog.fit <- function(Rmat, effects, Cmat, lb, ub, qp_pars){

  #----- Construct Cmat and bvec from lb and ub
  if (NROW(Cmat) > 0){
    # Get equality constraints
    iseq <- lb == ub
    meq <- sum(iseq)
    Amat <- Cmat[iseq,, drop = F]
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
  } else {
    # Null matrix
    Amat <- Cmat
    bvec <- NULL
    meq <- 0
    cmap <- integer(0)
  }

  #----- Fit QP

  # Normalise R matrix in case of big numbers
  sc <- norm(Rmat, "2")
  Rmat2 <- Rmat / sqrt(sc)
  effects2 <- effects / sqrt(sc)

  # Compute d vector
  dvec <- crossprod(effects2, Rmat2)

  # Fit QP
  res <- quadprog::solve.QP(solve(Rmat2), dvec, t(Amat), bvec, meq,
    factorized = TRUE)

  # Extract active constraints
  iact <- cmap[res$iact]

  # Return
  list(solution = res$solution, iterations = res$iterations[1],
    iact = sort(unique(iact)))
}
