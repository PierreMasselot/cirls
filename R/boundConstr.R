################################################################################
#
# Method to create a constraint matrix for constraint on the bound of smooths
# Works with various type of splines found in R.
#
################################################################################

#' Boundary constraints
#'
#' @description
#' Generic function to build a constraint matrix and bound vectors for constraining the boundary of [B-splines][boundConstr.bs()] and other nonlinear functions to a given value. It is designed to be used in the [constr][buildCmat()] formula interface. It allows methods for a wide range of regression terms.
#'
#' @param x An object representing a design matrix of predictor variables, typically basis functions. See Details for supported objects.
#' @param value A numerical indicating the boundary constraint value. Default to 0.
#' @param side Character indicating the side on which the constraint applies. One of `"right"` (the default), `"left"` or `"both"`.
#' @param thr For the factor method, an integer indicating a level from which the association should be equal to `value`. Cannot be used with `side = both`. This parameter has a slightly different interpretation for other methods.
#' @param sm A positive integer indicating the degree of smoothness for the constraint. In the default method indicates the number of coefficients to constrain to `value`, but can have a slightly different interpretation in other methods.
#' @param intercept For the default and factor method, a logical value indicating if the design matrix includes an intercept. In most cases, it will be automatically extracted from `x`, but this argument can be used to override it.
#' @param ... Additional parameters passed to or from other methods. Includes `value` and `side` for most methods.
#'
#' @details
#' Enforcing boundary constraints amounts to imposing an equality constraint on the `sm` first (if `side = left`), last (if `side = right`), or both (if `side = both`) coefficients of the basis matrix `x`. The equality constraint sets the bounds to `lb = ub = value`.
#'
#' ## Usage
#'
#' The recommended usage is to use this function through a call to `bound` on a term in the [constr][buildCmat()] formula interface. This method is then called internally to create the constraint matrix and bound vectors. However, `boundConstr` can also be called directly on a matrix-like object to manually build or inspect the constraint matrix.
#'
#' All methods internally rely on the default method for general matrices. Unless specified, all methods use the same parameters as the default one, which are passed through `...`. The only exception is `intercept` which is often inferred from the `x` object itself. In a typical usage in which `boundConstr` is called from the `constr` argument, `intercept` is automatically determined from the [glm][stats::glm()] formula.
#'
#' ## The `sm` parameter
#'
#' The `sm` parameter indicates the number of coefficients on the left/right that are constrained to be equal to `value` and can be interpreted as a smoothness degree for the boundary constraint. The default is to constraint only one coefficient. Note that `sm` has a slightly different interpretation when `x` represents [spline bases][boundConstr.bs()].
#'
#' ## Available methods
#'
#' In addition to the default method, `boundConstr` currently supports methods for several classes. The full list can also be consulted through `methods(boundConstr)`.
#'
#' ### Categorical variables
#'
#' * [factor][factor()]: for `factor` objects. Extract the [contrasts][stats::contrasts()] to define the constraint matrix. Here the `intercept` argument has the same interpretation as in the default method, i.e. if set to `TRUE` it means the `glm` model does not include an intercept externally to the factor. Note that, in this case, a simple dummy coding is done in R.
#' * [strata][boundConstr.strata()]: Indicator variables defining strata from the [dlnm][dlnm::dlnm()] package.
#'
#' ### Splines
#'
#' * [bs, ns, ps][boundConstr.bs()]: Specify boundary constraints on various flavors of B-splines.
#'
#' ### Distributed-lag linear and nonlinear models (DLNM)
#'
#' * [onebasis][boundConstr.onebasis()]: General method for basis functions generated in the package. Internally calls the relevant method for the specified basis.
#' * [crossbasis][boundConstr.crossbasis()]: Most exhaustive method to constrain DLNMs. See its dedicated [help page][boundConstr.crossbasis()].
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub`, respectively).
#'
#' @seealso [buildCmat][buildCmat()] detailing the `constr` interface.
#' @aliases bound
#'
#' @example inst/examples/ex_warming_bound_factor.R
#'
#' @order 1
#' @export
boundConstr <- function(x, ...) UseMethod("boundConstr")
