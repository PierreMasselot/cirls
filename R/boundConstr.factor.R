################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# factor method
#
################################################################################

#' @rdname boundConstr
#' @order 3
#' @export
boundConstr.factor <- function(x, value = 0, side = "right", thr = NULL, sm = 1,
  intercept = FALSE, ...){

  # Get the design matrix
  xmat <- stats::model.matrix(~ x)
  df <- NCOL(xmat)

  # Get levels to be constrained
  if (!is.null(thr)){
    sm <- switch(side,
      right = df - thr + 1,
      left = thr,
      stop("When 'thr' is provided, 'side' should be either 'right' or 'left'"))
  }

  # Get initial constraint matrix from default methods
  cm <- boundConstr.default(xmat, value = value, side = side, sm = sm,
    intercept = TRUE, ...)

  # Apply contrast if intercept is not included
  # NB: no contrast is applied in R when the model does not include an intercept
  if (!intercept & any(cm$Cmat[,1] != 0)) warning(
    "Need an intercept included to bound constrain on the left")
  if (isFALSE(intercept)){
    ctr <- stats::contrasts(x)
    cm$Cmat <- cm$Cmat %*% ctr
    cm$lb <- rep_len(cm$lb, NROW(cm$Cmat))
    cm$ub <- rep_len(cm$ub, NROW(cm$Cmat))
  }

  # Return
  cm
}
