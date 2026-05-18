################################################################################
#
# Reduce constraint matrix by removing constraints that can be removed safely
#
################################################################################

#' Reduce a set of linear constraints
#'
#' @description
#' Internal function checking whether a set of `Cmat`, `lb` and `ub` contains redundant or underlying equality constraints, and attempting to reduce the set accordingly. The objective is to obtain the smallest possible set of constraints without affecting the original problem.
#'
#' @param Cmat A constraint matrix.
#' @param lb,ub Bound vectors. If provided, should be consistent with `Cmat`. if not, default to a vector of 0 and `Inf`, respectively.
#' @param redundant Logical switch indicating whether to check for redundant constraints.
#' @param equality Logical switch indicating whether to check for underlying equality constraints.
#' @param warn Logical indicating if the user should be warned when constraint reduction happens.
#'
#' @details
#' This function is mostly for internal used, being silently called by [buildCmat][buildCmat()] and in some Constr functions used to create `Cmat`/`lb`/`ub`. Redundant or underlying equality constraints can result from putting together several sets of constraints together (for example non-decreasing and convex regression) and reducing the constraint matrix can be crucial to allow for inference. Therefore, if this check can be turned off when using `cirls` (see [cirls.control][cirls.control()]), it is not recommended to do so.
#'
#' ## Redundant constraint
#'
#' In a set of linear constraints \eqn{\mathbf{C}\beta \geq \mathbf{l}}, the constraint k is **redundant** if it can be removed without affecting the feasible region of the regression. Redundancy can be checked by finding the minimum possible value for the constraint \eqn{\mathbf{c}_k^T \beta} subject to all other constraints, i.e. by solving the linear optimisation problem
#'
#' \deqn{\begin{aligned}
#' \min \mathbf{c}_k^T \beta \\
#' s.t. \mathbf{C}_{-k}\beta \geq \mathbf{l}_{-k}
#' \end{aligned}}
#'
#' where \eqn{\mathbf{C}_{-k}} and \eqn{l_{-k}} are the constraint matrix and bound vector to which the kth row has been removed.
#'
#' If the minimum possible value is greater than or equal to \eqn{l_k}, then the corresponding constraint is redundant. This is easily extended for upper bounds, replacing the min by max in the problem above, and checking the maximum against the upper bound \eqn{u_k}.
#'
#' ## Underlying equality
#'
#' There is an underlying equality constraint when a constraint is forced to be equal to its lower bound \eqn{l_k} by other constraints. If the maximum possible value for \eqn{\mathbf{c}_k^T \beta} subject to all other constraints is equal to \eqn{l_k}, then it means \eqn{\mathbf{c}_k^T \beta} is forced to be exactly \eqn{l_k} and is therefore an underlying equality constraint. This is checked by changing the min to max in the problem above.
#'
#' ## Algorithm
#'
#' This function implements a simple algorithm to reduce the set of constraints. There is an initial check for row of `Cmat` that only contain zeros, as well as for constraints with infinite bounds only. Then it goes through every constraint with the following steps:
#'
#' 1. Find the minimum and maximum possible values of \eqn{\mathbf{c}_k^T \beta} using the optimisation problem above.
#' 2. If the minimum is higher than \eqn{l_k}, remove the bound, and same if the maximum is lower than \eqn{u_k}.
#' 3. If no bound remains, this constraint is discarded.
#' 4. If bounds remain, check for underlying equality constraint by comparing the minimum/maximum to \eqn{u_k/l_k}.
#'
#' This algorithm is similar to the "naive" algorithm of Caron et al. (1989).
#'
#' @note The result can depend on the ordering of `Cmat` because constraints are checked in the order provided. However, the results will not impact the feasible region in `cirls`, only what the matrix looks like.
#'
#' @returns A list with the following elements:
#' \item{Cmat/lb/ub}{The reduced set of constraints.}
#' \item{redundant}{Vector of indices indicating the redundant constraints that have been removed.}
#' \item{equality}{Vector of indices indicating the constraints that were part of underlying equality constraints.}
#'
#' @seealso [buildCmat][buildCmat()] for how constraint matrices are built in `cirls`.
#'
#' @references
#' Caron, R.J., McDonald, J.F., Ponic, C.M., 1989. A degenerate extreme point strategy for the classification of linear constraints as redundant or necessary. *J Optim Theory Appl* **62**, **225–237**. [DOI:10.1007/BF00941055](https://doi.org/10.1007/BF00941055)
#'
#' Telgen, J., 1983. Identifying Redundant Constraints and Implicit Equalities in Systems of Linear Constraints. *Management Science* **29**, **1209–1222**. [DOI:10.1287/mnsc.29.10.1209](https://doi.org/10.1287/mnsc.29.10.1209)
#'
#' @example inst/examples/ex_reduceCons.R
#'
#' @export
reduceCons <- function(Cmat, lb = NULL, ub = NULL,
  redundant = TRUE, equality = TRUE, warn = FALSE){

  # Dimensions
  m <- NROW(Cmat)
  p <- NCOL(Cmat)

  # Default values
  lb <- lb %||% rep(0, m)
  ub <- ub %||% rep(Inf, m)

  # Initialise an equality vector
  iseq <- rep(FALSE, m)

  # Start by checking if there are "zero" constraints which are useless
  # I use `all.equal` which takes a more sensible approach to equality to 0
  allzero <- apply(Cmat, 1, function(ck) isTRUE(all.equal(ck, rep(0, p))))
  lb[allzero] <- -Inf
  ub[allzero] <- Inf

  # Check if there "dummy" constraints with only infinite bounds and skip
  dummy <- lb == -Inf & ub == Inf
  inds <- which(!(dummy))

  # Only perform the following if either redundant or equality is on
  if (redundant | equality){

    # Avoid printing unwanted messages from linp
    sink(nullfile())

    # Go through all remaining constraints to perform "controlling" Lp
    for (k in inds){

      # Check if there is any non-redundant constraint left
      indup <- lb > -Inf | ub < Inf
      if (sum(indup[-k]) == 0) break

      # Build feasible region minus current constraint
      Cr <- rbind(Cmat[-k,, drop = F], -Cmat[-k,, drop = F])
      lr <- c(lb[-k], -ub[-k])
      keep <- lr < Inf & lr > -Inf

      #----- Check redundancy and implied equality

      # Min possible value of the cons
      lres <- limSolve::linp(Cost = Cmat[k,],
        G = Cr[keep,, drop = F], H = lr[keep], ispos = FALSE)

      # Max possible value
      ures <- limSolve::linp(Cost = -Cmat[k,],
        G = Cr[keep,, drop = F], H = lr[keep], ispos = FALSE)

      # Check redundancy
      # As for zeros above, consider a tolerance for small numerical errors
      if (!lres$IsError & redundant &
          lres$solutionNorm >= (lb[k] - sqrt(.Machine$double.eps))) lb[k] <- -Inf
      if (!ures$IsError & redundant &
          -ures$solutionNorm <= (ub[k] + sqrt(.Machine$double.eps))) ub[k] <- Inf

      # If both are redundant go to next constraint
      if (lb[k] == -Inf & ub[k] == Inf) next

      # Otherwise check for implied equality
      if (!lres$IsError & lres$solutionNorm == ub[k] & equality){
        iseq[k] <- TRUE
        lb[k] <- ub[k]
        next
      }
      if (!ures$IsError & -ures$solutionNorm == lb[k] & equality){
        iseq[k] <- TRUE
        ub[k] <- lb[k]
      }
    }

    # Remove the sink
    sink()
  }

  #----- Clean constraints

  # Exclude those that have been removed
  redund <- lb == -Inf & ub == Inf

  # Warn user
  rem <- sort(which(redund))
  if (warn & length(rem) > 0){
    warning(paste0("Redundant constraints removed: ",
      paste(rem, collapse = ", ")))
  }

  # Also warn about equality constraints
  eq <- sort(which(iseq))
  # if (warn & length(eq) > 0){
  #   warning(paste0("Underlying equality constraints: ",
  #     paste(eq, collapse = ", ")))
  # }

  # Return
  list(Cmat = Cmat[!redund,, drop = FALSE], lb = lb[!redund], ub = ub[!redund],
    redundant = rem, equality = eq)
}
