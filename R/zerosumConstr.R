################################################################################
#
# Method to create a zero-sum constraint matrix
# Useful in various cases, most notably for compositional regression
#
################################################################################

#' Zero-sum constraint matrix
#'
#' @description
#' Build constraint matrix and bounds for coefficients summing to zero.
#'
#' @param ... Variables to be included in the constraint.
#' @param group If set to TRUE, the constraint is build independently for each variable in `...`
#'
#' @export
zerosumConstr <- function(..., group = FALSE){

  # Extract
  varlist <- list(...)

  # Extract number of columns in each
  nv <- length(varlist)
  ncs <- sapply(varlist, NCOL)
  nctot <- sum(ncs)

  # Create constraint matrix depending on groups
  if (isFALSE(group)){
    Cmat <- t(rep(1, nctot))
    attr(Cmat, "lb") <- 0
    attr(Cmat, "ub") <- 0
  } else {
    Cmat <- matrix(0, nv, nctot)
    Cmat[cbind(rep(seq_len(nv), ncs), seq_len(nctot))] <- 1
    attr(Cmat, "lb") <- rep(0, nv)
    attr(Cmat, "ub") <- rep(0, nv)
  }

  # Return
  Cmat
}
