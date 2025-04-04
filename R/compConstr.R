################################################################################
#
# Method to create a constraint matrix for compositional regression
#
################################################################################

compConstr <- function(x){

  # Check data
  iscomp <- all.equal(rowSums(exp(x)), rep(1, NROW(x)))
  if (!iscomp) warning("The data provided in x might not be log-compositional.")

  # Extract info on data
  p <- ncol(x)

  # Create constraint matrix and return with bounds
  list(Cmat = t(rep(1, p)), lb = 0, ub = 0)
}
