################################################################################
#
# Create full derivative matrix for B-splines
# Not exported
#
################################################################################

#' Spline derivative matrix
#'
#' @description
#' Computes a derivative matrix for B-splines that can then be used for shape-constraints. It is internally called by [shapeConstr][shapeConstr()] and should not be used directly.
#'
#' @param d Non-negative integer giving the order of derivation. Should be between 0 and `ord - 2`.
#' @param s Sign of the derivative.
#' @param knots Vector of ordered knots from the spline bases.
#' @param ord Non-negative integer giving the order of the spline.
#'
#' @details Does the heavy lifting in [shapeConstr][shapeConstr()] to create a constraint matrix for shape-constrained B-splines. Only useful for advanced users to create constraint matrices without passing an object to one of the [shapeConstr][shapeConstr()] methods.
#'
#' @note
#' `dmat` doesn't perform any checks of the parameters so use carefully. In normal usage, checks are done by [shapeConstr][shapeConstr()] methods.
#'
#' @returns A matrix of weighted differences that can be used to constrain B-spline bases.
#'
#' @examples
#' # A second derivative matrix for cubic B-Splines with regularly spaced knots
#' # Can be used to enforce convexity
#' cirls:::dmat(2, 1, 1:15, 4)
#'
dmat <- function(d, s, knots, ord){

  # Create diff matrices for each derivative order (including "zero" deriv)
  alldm <- lapply(rev(seq_len(d)), dm, knots = knots, ord = ord)
  alldm <- c(alldm, list(diag(length(knots) - ord)))

  # Matrix multiply everything
  s * Reduce("%*%", alldm)
}
