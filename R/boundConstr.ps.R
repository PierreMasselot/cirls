################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# dlnm:::ps method
#
################################################################################

# By default deg is set to the degree of the spline, for smoother convergence

#' @rdname boundConstr.bs
#' @order 3
#' @export
boundConstr.ps <- function(x, ...)
{

  # Check if a degree is provided
  degree <- attr(x, "degree")
  dots <- list(...)
  if ("sm" %in% names(dots)) if(dots$sm < degree) stop(
    "'sm' needs to be at least equal to the degree of the spline")

  # Get constraint matrix
  boundConstr.bs(x, ...)
}
