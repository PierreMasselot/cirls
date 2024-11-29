################################################################################
#
# Create full derivative matrix for B-splines
#
################################################################################

dmat <- function(l, knots, ord){
  # Create diff matrices for each derivative order
  alldm <- lapply(l:1, dm, knots = knots, ord = ord)

  # Matrix multiply everything
  Reduce("%*%", alldm)
}
