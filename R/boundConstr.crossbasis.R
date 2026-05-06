################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# crossbasis method
#
################################################################################

#' @rdname shapeConstr.crossbasis
#' @export
boundConstr.crossbasis <- function(x, dim = NULL, overall = FALSE, slice = NULL,
  ...){

  # Call cbConstr
  cmlist <- cbConstr(x, constr = "bound", pars = list(...),
    dim = dim, slice = slice, overall = overall)

  # Change bound and return
  cmlist$ub <- cmlist$lb
  cmlist
}
