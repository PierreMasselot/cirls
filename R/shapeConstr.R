################################################################################
#
# Method to create a constraint matrix for shape-constrained splines
# Works with various type of splines found in R.
#
################################################################################

#' Specify shape constraints
#'
#' @description
#' Builds a constraint matrix and bound coefficients for shape-constraints on a set of coefficients. This is a generic function designed to be used in the [constr][buildCmat()] interface. It allows methods for a wide range of regression terms (see details).
#'
#' @param x An object representing a design matrix of predictor variables, typically basis functions. See details for supported objects.
#' @param shape A character vector indicating one or several shape-constraints. See details for supported shapes.
#' @param intercept For the default method, a logical value indicating if the design matrix includes an intercept. In most cases, it will be automatically extracted from `x` but this can be used to override it.
#' @param ... Additional parameters passed to or from other methods.
#'
#' @details
#' This function is used to specify shape constraints on terms in a `cirls` model. Shapes can refer to the relation between the coefficients of several variables in a matrix (for instance dummy variables to impose shapes on the levels of a factor), or to the shape of a smooth term (such as splines). note that this function can also be used to easily specify non-negativity or non-positivity constraints. See below for the list of implemented shapes and the examples section for several use cases.
#'
#' ## Usage
#'
#' The recommended usage is to use this function through a call to `shape` on a term in the [constr][buildCmat()] interface. This method is then called internally to create the constraint matrix and bound vectors. However, `shapeConstr` can also be called directly on a matrix-like object to manually build or inspect the constraint matrix.
#'
#' The parameters necessary to build the constraint matrix (e.g. `knots` and `ord` for splines) are typically extracted from the `x` object. This is also true for the `intercept` for most of the object, except for the default method for which it can be useful to explicitly provide it. In a typical usage in which `shapeConstr` is called from the `constr` argument, `intercept` is automatically determined from the [glm][stats::glm()] formula.
#'
#' ## Allowed shapes
#'
#' The `shape` argument allows to define a specific shape for the association between the expanded term in `x` and the response of the regression model. This shape can describe the relation between coefficients for the default method, or the shape of the smooth term for spline bases. At the moment, six different shapes are supported, with up to three allowed simultaneously (one from each category):
#'
#' * `"pos"` or `"neg"`: Positive/Negative. Applies to the full association.
#' * `"inc"` or `"dec"`: Monotonically Increasing/Decreasing.
#' * `"cvx"` or `"ccv"`: Convex/Concave.
#'
#' ## Available methods
#'
#' In addition to the default method, `shapeConstr` currently supports several objects, creating an appropriate shape-constraint matrix depending on the object. The full list can be obtained by `methods(shapeConstr)`.
#'
#' ### General
#' * [factor()]: for categorical variables. Extract the [contrasts][stats::contrasts()] to define the constraint matrix. Here the `intercept` argument has the same interpretation as in the default method, i.e. if set to `TRUE` it means the `glm` model doesn't include an intercept externally to the factor. Note that, in this case, a simple dummy coding is done in R.
#'
#' ### From the [splines][splines::splines] package
#'
#' * [bs][splines::bs()]: B-splines.
#' * [ns][splines::ns()]: Natural splines.
#'
#' ### From the [dlnm][dlnm::dlnm] package
#'
#' * [onebasis][dlnm::onebasis()]: General method for basis functions generated in the package.
#' * [ps][dlnm::ps()]: Penalised splines (P-Splines).
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub` respectively).
#'
#' @references
#' Zhou, S. & Wolfe, D. A., 2000. On derivative estimation in spline regression. *Statistica Sinica* **10**, **93–108**.
#'
#' @seealso [buildCmat][buildCmat()] detailing the `constr` interface.
#'
#' @example inst/examples/ex_london_nonneg.R
#' @example inst/examples/ex_warming_factor.R
#' @example inst/examples/ex_warming_splines.R
#'
#' @export
shapeConstr <- function(x, shape, ...) UseMethod("shapeConstr")
