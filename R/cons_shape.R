################################################################################
#
# Method to create a constraint matrix for shape-constrained splines
# Works with various type of splines found in R.
#
################################################################################

#' Constraint matrices for shape-constrained splines
#'
#' @export
cons_shape <- function(x, ...) UseMethod("cons_shape")
