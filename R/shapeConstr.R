################################################################################
#
# Method to create a constraint matrix for shape-constrained splines
# Works with various type of splines found in R.
#
################################################################################

#' Constraint matrices for shape-constrained splines
#'
#' @export
shapeConstr <- function(x, shape, ...) UseMethod("shapeConstr")
