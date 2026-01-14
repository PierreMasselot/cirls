################################################################################
#
# Method to create a zero-sum constraint matrix
# Useful in various cases, most notably for compositional regression
#
################################################################################

#' Specify Zero-sum constraints
#'
#' @description
#' Builds a constraint matrix and bound vectors for zero-sum constraints on a set of coefficients. This is a generic function designed to be used in the [constr][buildCmat()] interface.
#'
#' @param ... Terms to be included in the constraint.
#' @param group If set to TRUE, the constraint is build independently for each variable in `...`
#'
#' @details
#' This function is used to force one or several set(s) of coefficients to sum to zero. The main use is for when (log) compositions are used as predictors in the model, but it can also be used more generally for variables that are relative to a reference.
#'
#' ## Usage
#'
#' The recommended usage is to use this function through a call to `zerosum` on a term in the [constr][buildCmat()] interface. This method is then called internally to create the constraint matrix and bound vectors. However, `shapeConstr` can also be called directly on a matrix-like object to manually build or inspect the constraint matrix.
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub` respectively).
#'
#' @references
#' Altenbuchinger, M., Rehberg, T., Zacharias, H.U., Stämmler, F., Dettmer, K., Weber, D., Hiergeist, A., Gessner, A., Holler, E., Oefner, P.J., Spang, R., 2017. Reference point insensitive molecular data analysis. *Bioinformatics* **33**, **219–226**. [DOI:10.1093/bioinformatics/btw598](https://doi.org/10.1093/bioinformatics/btw598)
#'
#' Aitchison, J., Bacon-Shone, J., 1984. Log contrast models for experiments with mixtures. Biometrika 71, 323–330. [DOI:10.1093/biomet/71.2.323](https://doi.org/10.1093/biomet/71.2.323)
#'
#' @seealso [buildCmat][buildCmat()] detailing the `constr` interface.
#'
#' @example inst/examples/ex_fgl_coda.R
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
  } else {
    Cmat <- matrix(0, nv, nctot)
    Cmat[cbind(rep(seq_len(nv), ncs), seq_len(nctot))] <- 1
  }

  # Bounds
  lb <- rep(0, NROW(Cmat))
  ub <- rep(0, NROW(Cmat))

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
