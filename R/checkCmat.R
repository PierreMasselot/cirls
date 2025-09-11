################################################################################
#
# Check if constraint matrix is reducible
# Returns constraints that can be removed safely
#
# Loop over constraints:
#   1. Check if the constraint can be expressed as positive linear combination
#     of other constraints
#     - This is done by the `coneB` routine (a non-negative linear least-squares)
#   2. If so, check if the redundancy represents an equality constraint
#     - This is done by checking is the current constraint is the opposite of a
#       positive linear combination of other constraints
#
# NB: coneB solves NNLS problems
#
################################################################################

#' Check constraint matrix irreducibility
#'
#' @description
#' Checks a constraint matrix does not contains redundant rows
#'
#' @param Cmat A constraint matrix as passed to [cirls.fit()]
#'
#' @details
#' The user typically doesn't need to use `checkCmat` as it is internally called by [cirls.control()]. However, it might be useful to undertsand if `Cmat` can be reduced for inference purpose. See the note in [confint.cirls()].
#'
#' A constraint matrix is irreducible if no row can be expressed as a *positive* linear combination of the other rows. When it happens, it means the constraint is actually implicitly included in other constraints in the matrix and can be dropped. Note that this a less restrictive condition than the constraint matrix having full row rank (see some examples).
#'
#' The function starts by checking if some constraints are redundant and, if so, checks if they underline equality constraints. In the latter case, the constraint matrix can be reduced by expressing these constraints as a single equality constraint with identical lower and upper bounds (see [cirls.fit()]).
#'
#' The function also checks whether there are "zero constraints" i.e. constraints with only zeros in `Cmat` in which case they will be labelled as redundant.
#'
#' @returns A list with three elements:
#' \item{redundant}{Logical vector of indicating redundant constraints}
#' \item{equality}{Logical vector indicating which constraints are part of an underlying equality constraint}
#'
#' @seealso [confint.cirls()]
#'
#' @references
#' Meyer, M.C., 1999. An extension of the mixed primal–dual bases algorithm to the case of more constraints than dimensions. *Journal of Statistical Planning and Inference* **81**, 13–31. \doi{10.1016/S0378-3758(99)00025-7}
#' @example man/examples/checkCmat_ex.R
#'
#' @export
checkCmat <- function(Cmat){
  # Check if there are "zero" constraints which are redundant
  # I use `all.equal` which takes a more sensible approach to equality to 0
  redundant <- apply(Cmat, 1,
    function(x) isTRUE(all.equal(x, rep(0, ncol(Cmat)))))
  # Prepare Cmat
  tCmat <- t(Cmat)
  equality <- rep(FALSE, ncol(tCmat))
  for (i in which(!redundant)){
    # Break the loop if there is only one useful constraint left
    if (sum(!redundant) < 2) break
    # Check redundancy
    y <- tCmat[, i]
    x <- tCmat[, -c(i, which(redundant)), drop = F]
    # res <- coneproj::coneB(y, x)$yhat # Returns error for some problems
    fit <- limSolve::nnls(x, y)
    res <- x %*% fit$X
    redundant[i] <- isTRUE(all.equal(y, drop(res)))
    if (!redundant[i]){
      # Check underlying equality constraint
      fiteq <- limSolve::nnls(x, -y)
      reseq <- x %*% fiteq$X
      # reseq <- coneproj::coneB(-y, x)$yhat # Returns error for some problems
      equality[i] <- isTRUE(all.equal(-y, drop(reseq)))
    }
  }
  # Return indices
  list(redundant = redundant, equality = equality)
}

