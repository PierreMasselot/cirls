################################################################################
#
# Shape constraint matrix method:
# Default method for a general matrix
#
################################################################################

#' @rdname shapeConstr
#' @order 2
#' @export
shapeConstr.default <- function(x, shape, range = NULL, intercept = FALSE, ...) {

  #----- Parameters

  # Matrix dimension and check
  df <- NCOL(x) + !intercept
  cpars <- chkshp(shape, df)

  # Check range parameter
  rng <- chkrng(range, 1:df, msg = FALSE)

  #----- Create matrix

  # Create constraint matrices
  Cmat <- lapply(cpars, function(cp){
    m <- diag(df)
    if (cp[1] > 0) m <- diff(m, diff = cp[1])
    cp[2] * m
  })

  # Put together
  Cmat <- do.call(rbind, Cmat)

  # Remove intercept if not included
  if (!intercept) Cmat <- Cmat[,-1, drop = F]

  # Keep only desired range
  rind <- (1:ncol(Cmat)) >= rng[1] & (1:ncol(Cmat)) <= rng[2]
  keep <- apply(Cmat, 1, function(r) all(r[!rind] == 0))
  Cmat <- Cmat[keep,, drop = FALSE]

  # Reduce if necessary
  Cmat <- reduceCons(Cmat, warn = FALSE)$Cmat

  # Add bounds attributes
  clist <- list(Cmat = Cmat, lb = rep(0, NROW(Cmat)), ub = rep(Inf, NROW(Cmat)))

  # Return
  clist
}
