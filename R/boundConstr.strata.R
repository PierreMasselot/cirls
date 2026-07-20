################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# strata method
#
################################################################################

#' @rdname boundConstr.onebasis
#' @order 2
#' @export
boundConstr.strata <- function(x, value = 0, side = "right", thr = NULL, sm = 1,
  ...){

  # Get info
  int <- attr(x, "intercept")
  ref <- attr(x, "ref")
  breaks <- attr(x, "breaks")
  ncat <- length(breaks) + 1

  # Get the number of levels to constrain
  if (!is.null(thr)){
    sm <- switch(side,
      right = sum((breaks - thr) > sqrt(.Machine$double.eps)) + 1,
      left = sum((breaks - thr) < -sqrt(.Machine$double.eps)) + 1,
      stop("When 'thr' is provided, 'side' should be either 'right' or 'left'"))
  }

  # Get initial constraint matrix from default methods
  cm <- boundConstr.default(diag(ncat), value = value, side = side, sm = sm,
    intercept = TRUE, ...)

  # Adjust for reference: same as in dlnm:::strata
  if (!int & any(cm$Cmat[,1] != 0)) warning(
    "Need an intercept included to bound constrain on the left")
  if (ref > 0){
    cm$Cmat <- cm$Cmat[, -ref, drop = FALSE]
    if (int) cm$Cmat <- cbind(1, cm$Cmat)
  }

  # Return
  cm
}
