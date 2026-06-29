################################################################################
#
# Function to fit model by osqp
#
################################################################################

osqp.fit <- function(Rmat, effects, Cmat, lb, ub, qp_pars){

  # P and q
  Pmat <- crossprod(Rmat)
  qvec <- -crossprod(effects, Rmat)

  # Fit
  res <- osqp::solve_osqp(P = Pmat, q = qvec, A = Cmat, l = lb, u = ub,
    pars = qp_pars)

  # Return
  list(solution = res$x, iterations = res$info$iter,
    iact = which(res$y != 0))
}
