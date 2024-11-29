################################################################################
#
# Check parameters for the cons_shape method
#
################################################################################

chkcpars <- function(diff, sign, intercept, ord = Inf) {
  
  # Recycle vectors
  ncons <- max(lengths(list(diff, sign)))
  diff <- rep_len(diff, ncons)
  sign <- rep_len(sign, ncons)
  
  # Check diff parameter(s)
  diff <- unique(as.integer(diff))
  if (any(diff < 0 | diff >= ord)){
    stop("'diff' must be between 0 and the order of the spline")
  }
  
  # Check sign parameter(s)
  sign <- sign(as.numeric(sign))
  if (any(sign == 0)){
    stop("'sign' must be different than 0")
  }
  
  # Check intercept
  intercept <- as.logical(intercept)
  if (is.na(intercept)) stop("'intercept' must be TRUE or FALSE")
  
  # Return
  list(diff = diff, sign = sign, intercept = intercept)
}
