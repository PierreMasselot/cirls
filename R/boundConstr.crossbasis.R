################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# crossbasis method
#
################################################################################

#' @rdname shapeConstr.crossbasis
#' @export
boundConstr.crossbasis <- function(x, dim = "var", overall = FALSE, odrng = NULL,
  ...){

  # Call cbConstr
  cmlist <- cbConstr(x, constr = "bound", pars = list(...),
    dim = dim, odrng = odrng, overall = overall)

  # Change bound and return
  # cmlist$ub <- cmlist$lb
  cmlist
}
