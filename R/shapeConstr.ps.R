################################################################################
#
# Shape constraint matrix method:
# ps method
#
################################################################################

#' @rdname shapeConstr
#' @export
shapeConstr.ps <- function(x, shape, ...){

  # Same as B-Splines
  shapeConstr.bs(x, shape, ...)
}
