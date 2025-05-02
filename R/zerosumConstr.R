################################################################################
#
# Method to create a zero-sum constraint matrix
# Useful in various cases, most notably for compositional regression
#
################################################################################

#' @export
zerosumConstr <- function(x){

  # Create constraint matrix and return with bounds
  list(Cmat = t(rep(1, ncol(x))), lb = 0, ub = 0)
}
