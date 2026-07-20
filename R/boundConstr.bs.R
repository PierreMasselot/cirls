################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# bs method
#
################################################################################

#' Boundary constraints for spline bases
#'
#' @description
#' Method for constraining a smooth functions estimated using B-Spline bases to reach a given value at its boundary.
#'
#' @param x A B-spline basis object.
#' @param value A numerical indicating the boundary constraint value. Default to 0.
#' @param side Character indicating the side on which the constraint applies. One of `"right"` (the default), `"left"` or `"both"`.
#' @param thr A numeric value indicating a threshold from which the smooth should be equal to `value`. Default to the boundary knots of the spline basis. Cannot be used with `side = both`. See details.
#' @param sm A positive integer controlling the smoothness of the constraint. Ignored when `thr` is provided. See details.
#' @param ... Additional parameters passed to or from other methods.
#'
#' @details
#' See the [generic][boundConstr()] method for an overview of bound constraints.
#'
#' ## Boundary constraints with B-splines
#'
#' By construction, B-spline basis functions sum to one at any point. Therefore, boundary constraints can essentially be enforced through equality constraints on the first (if `side = left`) or last (if `side = right`) few coefficients of the B-spline basis, with small variations for [natural][splines::ns()] or [penalised][dlnm::ps()] splines.
#'
#' ## The `thr` and `sm` arguments
#'
#' The `sm` argument can control how smoothly the resulting function reaches the constraint `value`. If only one basis is non-null at the boundary and need its coefficient to be constrained for the boundary constraint, the number of bases that reach zero exactly at the boundary is equal to the degree of the spline. Forcing their coefficients to equal `value` results in a smoother convergence at the boundary. Note that increasing `sm` to be higher than the degree of the spline will propagate the equality constraint towards interior knots.
#'
#' Alternatively, the user can provide a threshold through `thr` above which (if `side = right`) or below which (if `side = left`) the resulting function is guaranteed to be equal to `value`. This allows propagating the boundary constraint towards the interior of the range. Note that where the constraint starts depends on the knots, as using `thr` will actually enforce the constraint from the maximum (minimum) interior knot such as `thr` is included in the constrained range. This does not allow `side = both` so if thresholds are needed on both sides of the range, it needs to be done with two calls to `boundConstr`.
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub`, respectively).
#'
#' @seealso The generic [boundConstr][boundConstr()] functions.
#'
#' @example inst/examples/ex_warming_bound_bs.R
#'
#' @order 1
#' @export
boundConstr.bs <- function(x, value = 0, side = "right", thr = NULL,
  sm = NULL, ...)
{

  # Get info from basis
  degree <- attr(x, "degree")
  df <- ncol(x)
  int <- attr(x, "intercept")
  ik <- attr(x, "knots")
  bk <- attr(x, "Boundary.knots")
  kn <- sort(c(rep(bk, degree + 1), ik))

  #----- Determine the number of bases to constraints

  # Sense check thr
  if (isTRUE(thr < kn[1] & thr > kn[length(kn)])) warning(
    paste0("'thr' should be",
    " between bounday knots. Constraining at the boundary only."))

  # If provided thr takes precedence
  if (!is.null(thr)){
    sm <- switch(side,
      right = sum((kn - thr) > sqrt(.Machine$double.eps)),
      left = sum((kn - thr) < -sqrt(.Machine$double.eps)),
      stop("When 'thr' is provided, 'side' should be either 'right' or 'left'"))
  } else {
    # If not provided then we use smoothing degree
    # Set by default at the degree of spline for smoother convergence
    sm <- sm %||% degree
    if (sm > degree) warning(paste0("Setting 'sm' greater than the degree",
      " of the spline will also constrain towards interior knots"))
  }

  # Get constraint matrix
  boundConstr.default(x, value = value, side = side, sm = sm, intercept = int,
    ...)
}
