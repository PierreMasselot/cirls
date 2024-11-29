################################################################################
#
# Shape constraint matrix method:
# ps method
#
################################################################################

#' @export
cons_shape.ps <- function(x, diff = 0, sign = 1, ...){

  # Same as B-Splines
  cons_shape.bs(x, diff, sign)
}
