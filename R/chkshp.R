################################################################################
#
# Check parameters for the shapeConstr method
#
################################################################################

# chkshp <- function(shape, ord = Inf) {
#
#   # Check constraints: first element is diff, second is sign
#   conslist <- list("pos" = c(0, 1), "neg" = c(0, -1), "inc" = c(1, 1),
#     "dec" = c(1, -1), "cvx" = c(2, 1), "ccv" = c(2, -1))
#   shape <- match.arg(shape, names(conslist), several.ok = TRUE)
#
#   # Extract differences and signs
#   cpars <- conslist[[shape]]
#
#   # Check that difference works with order
#   if (cpars[1] >= ord) stop(
#     "Constrained difference must be lower than the order of the spline")
#
#   # Return
#   cpars
# }


chkshp <- function(shape, ord = Inf) {

  # Check constraints: first element is diff, second is sign
  conslist <- list("pos" = c(0, 1), "neg" = c(0, -1), "inc" = c(1, 1),
    "dec" = c(1, -1), "cvx" = c(2, 1), "ccv" = c(2, -1))
  shape <- match.arg(shape, names(conslist), several.ok = TRUE)

  # Extract differences and signs
  cpars <- lapply(shape, function(x) conslist[[x]])

  # Check difference
  dvec <- sapply(cpars, "[", 1)
  wr <- dvec >= ord
  if (any(wr)) stop("Constrained difference must be lower than the order of the spline")

  # Check non-redundant shapes
  dup <- duplicated(dvec)
  if (any(dup)){
    cpars <- cpars[!dup]
    warning(sprintf("Non-compatible shapes provided, only [%s] kept",
      paste(shape[!dup], collapse = ", ")))
  }

  # Return
  cpars
}
