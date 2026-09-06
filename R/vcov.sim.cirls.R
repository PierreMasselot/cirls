################################################################################
#
#    Variance-covariance matrix method for sim.cirls objects
#
################################################################################

#' @rdname simulCoef
#' @order 5
#' @export
vcov.sim.cirls <- function(object, ...){

  # Compute the variance-covariance matrix directly
  stats::var(object)
}
