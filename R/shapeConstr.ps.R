################################################################################
#
# Shape constraint matrix method:
# ps method
#
################################################################################

#' @export
shapeConstr.ps <- function(x, shape, ...){

  # Same as B-Splines
  shapeConstr.bs(x, shape, ...)
}
