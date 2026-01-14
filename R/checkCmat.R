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

#' Check and reduce constraint matrix
#'
#' @description
#' Checks whether some constraints are redundant and can be removed, and whether some constraints can be reduced into an equality constraint. If requested, redundant constraints are removed, but not equality ones.
#'
#' @param Cmat A constraint matrix.
#' @param lb,ub Bound vectors. If provided, should be consistent with `Cmat`.
#' @param reduce If TRUE, returns reduced `Cmat`, `lb` and `ub` with redundant constraints removed. Otherwise, they are returned as provided.
#' @param warn Whether to warn the user when redundant and equality constraints are found.
#'
#' @details
#' The user typically doesn't need to use `checkCmat` as it is internally called by [buildCmat][buildCmat()] and in some Constr functions creating `Cmat`/`lb`/`ub` . However, it might be useful in the case of inconsistent constraints or underlying equality constraints.
#'
#' ## Irreducibility
#'
#' `Cmat` is irreducible if there is no *redundant* constraint, and no *underlying equality* constraint.
#'
#' A row in `Cmat` is *redundant* if it can be expressed as a *positive* linear combination of other rows. When it happens, it means the corresponding constraint is actually implicitly included in other constraints, and can be dropped without affecting the problem. Rows that only contain zeros are considered redundant.
#'
#' When it exists a linear combination of rows in `Cmat` that result in the null vector, it means that there is an *underlying equality* constraint. In that case, it means the corresponding rows can be reduced to a single equality constraint. Note that underlying constraints are left as is in the returned `Cmat` even when `reduce = TRUE`, and it is left to the user to figure out whether to reduce it.
#'
#' @note `checkCmat` works when only `Cmat` is provided. `lb`/`ub` can be provided to conveniently be reduced in the case of redundant constraints, but they should be consistent with `Cmat`.
#'
#' @returns A list with the following elements:
#' \item{redundant}{Logical vector indicating redundant constraints in the provided `Cmat`.}
#' \item{equality}{An integer vector indicating the groups of underlying equality constraints in the provided `Cmat`. Zero values indicate that the row is not part of any underlying equality constraint and higher values indicate which equality constraint it is part of.}
#' \item{Cmat/lb/ub}{The reduced constraints when `reduce = TRUE` or the provided constraints otherwise.}
#'
#' @seealso [buildCmat][buildCmat()]
#'
#' @references
#' Meyer, M.C., 1999. An extension of the mixed primal–dual bases algorithm to the case of more constraints than dimensions. *Journal of Statistical Planning and Inference* **81**, 13–31. [DOI:10.1016/S0378-3758(99)00025-7](https://doi.org/10.1016/S0378-3758(99)00025-7)
#'
#' @example inst/examples/ex_checkCmat.R
#'
#' @export
checkCmat <- function(Cmat, lb = NULL, ub = NULL, reduce = TRUE, warn = TRUE){
  # Check if there are "zero" constraints which are useless
  # I use `all.equal` which takes a more sensible approach to equality to 0
  redundant <- apply(Cmat, 1,
    function(x) isTRUE(all.equal(x, rep(0, ncol(Cmat)))))
  # Prepare Cmat
  tCmat <- t(Cmat)
  equality <- rep(0, ncol(tCmat))
  eq_cnt <- 0
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

    # Check equality only if not redundant and not already detected
    if (!redundant[i] & equality[i] == 0){
      # Check underlying equality constraint: origin is a linear combination
      fiteq <- limSolve::nnls(x, -y)
      reseq <- x %*% fiteq$X
      # reseq <- coneproj::coneB(-y, x)$yhat # Returns error for some problems

      # If there is an equality constraint, check which variables involved
      if (isTRUE(all.equal(-y, drop(reseq)))){
        eq_cnt <- eq_cnt + 1
        eqi <- c(i,
          setdiff(which(!redundant), i)[(fiteq$X - 0) >
              sqrt(.Machine$double.eps)])
        equality[eqi] <- eq_cnt
      }
    }
  }

  # Reduce constraints if requested
  if (any(redundant) && reduce){
    Cmat <- Cmat[!redundant,,drop = F]
    lb <- lb[!redundant]
    ub <- ub[!redundant]
  }

  # Warn user if requested
  if (warn){
    if (any(redundant)) warning(paste0("Redundant constraint found: ",
        paste(which(redundant), collapse = ", ")))
    if (eq_cnt > 0) warning(paste0(eq_cnt,
      " underlying equality constraint found"))
  }

  # Return indices
  list(redundant = redundant, equality = equality,
    Cmat = Cmat, lb = lb, ub = ub)
}

