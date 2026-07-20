################################################################################
#
# Bound constraint
# default method
#
################################################################################

#' @rdname boundConstr
#' @order 2
#' @export
boundConstr.default <- function(x, value = 0, side = c("right", "left", "both"),
  sm = 1, intercept = FALSE, ...)
{

  # Extract info
  df <- ncol(x)

  # Check side parameter
  side <- match.arg(side)

  # Check the number of constraints parameters
  if (sm > df | sm < 1) stop(
      "'sm' should be an integer between 1 and the number of bases")

  # Check intercept
  if (!intercept & side != "right") warning(
    "Need an intercept included to bound constrain on the left")

  # Create Cmat depending on side
  cr <- cl <- NULL
  if (side != "right"){
    ncl <- sm - 1 + intercept
    cl <- cbind(diag(ncl), matrix(0, ncl, df - ncl))
  }
  if (side != "left"){
    cr <- cbind(matrix(0, sm, df - sm), diag(sm))
  }
  Cmat <- rbind(cr, cl)

  # Return with bounds
  list(Cmat = Cmat, lb = rep(value, NROW(Cmat)), ub = rep(value, NROW(Cmat)))
}
