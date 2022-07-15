################################################################################
#
# Function to fit model by osqp
#
################################################################################

osqp.fit <- function(Dmat, dvec, Cmat, lb, ub, qp_pars){
  
  # Fit
  res <- osqp::solve_osqp(P = Dmat, q = -dvec, A = Cmat, l = lb, u = ub, 
    pars = qp_pars)
  
  # Return
  list(solution = res$x, iterations = res$info$iter, 
    iact = which(res$y != 0))
}