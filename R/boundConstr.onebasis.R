################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# onebasis method
#
################################################################################

#' Boundary constraints for onebasis objects
#'
#' @description
#' Method to set up boundary constraints for basis objects constructed through [onebasis][dlnm::onebasis()] and other specifications from the [dlnm][dlnm::dlnm()] package.
#'
#' @param x A `onebasis` object or any other objects from the `dlnm` package.
#' @param value A numerical indicating the boundary constraint value. Default to 0.
#' @param side Character indicating the side on which the constraint applies. One of `"right"` (the default), `"left"` or `"both"`.
#' @param thr An integer indicating a threshold from which the association should be equal to `value`. Using `thr` overrides `sm`. Cannot be used with `side = both`.
#' @param sm A positive integer indicating the degree of smoothness for the constraint. Typically indicates the number of coefficients to constrain to `value`.
#' @param ... Additional parameters passed to or from other methods.
#'
#' @details
#' The `onebasis` method simply check there is an available [boundConstr][boundConstr()] method available for the underlying basis type (for instance [strata][boundConstr.strata()] or [P-splines][boundConstr.ps()]) and dispatches `x` to the appropriate method. The ellipsis `...` is used to pass the necessary arguments to the method, including `value` and `side`. See the help for the specific methods.
#'
#' The `strata` method allows to specify boundary constraints for [strata][dlnm::strata()] objects, which are factors expanded through dummy parametrisation. In this context, boundary constraints essentially add an equality constraint to `value` on the `sm` extreme coefficients (depending on `side`). A numeric threshold `thr` can be specified instead of `sm`, which will add the equality constraint from the strata that contains `thr`. This allows to propagate boundary constraints towards the interior of the range.
#'
#' @return A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub`, respectively).
#'
#' @seealso The generic [boundConstr][boundConstr()] functions. The [crossbasis][boundConstr.crossbasis()] is used to specify the constraint on distributed lag (non-linear) models.
#'
#' @example inst/examples/ex_warming_bound_onebasis.R
#'
#' @order 1
#' @export
boundConstr.onebasis <- function(x, ...){

  # Extract the right method
  fun <- attr(x, "fun")
  met <- paste0("boundConstr.", fun)
  if (!met %in% utils::methods("boundConstr")) {
    warning(paste0("No existing 'boundConstr' method for '", fun,
      "' functions. Using default method."))
    met <- "boundConstr.default"
  }

  # Call the right method
  pars <- utils::modifyList(list(x = x), list(...))
  do.call(met, pars)
}
