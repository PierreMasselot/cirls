################################################################################
#
# Method to create a constraint matrix for shape-constrained splines
# Works with various type of splines found in R.
#
################################################################################

#' Create shape constraints
#'
#' @description
#' Creates a constraint matrix to shape-constrain a set of coefficients. Mainly intended for splines but can constrain various bases or set of variables. Will typically be called from within [cirls.fit][cirls.fit()] but can be used to generate constraint matrices.
#'
#' @param x An object representing a design matrix of predictor variables, typically basis functions. See details for supported objects.
#' @param shape A character vector indicating one or several shape-constraints. See details for supported shapes.
#' @param intercept For the default method, a logical value indicating if the design matrix includes an intercept. In most cases will be automatically extracted from `x`.
#' @param ... Additional parameters passed to or from other methods.
#'
#' @details
#' The recommended usage is to directly specify the shape constraint through the `shape` argument in the call to [glm][stats::glm()] with [cirls.fit][cirls.fit()]. This method is then called internally to create the constraint matrix. However, `shapeConstr` can nonetheless be called directly to manually build or inspect the constraint matrix for a given shape and design matrix.
#'
#' The parameters necessary to build the constraint matrix (e.g. `knots` and `ord` for splines) are typically extracted from the `x` object. This is also true for the `intercept` for most of the object, except for the default method for which it can be useful to explicitly provide it. In a typical usage in which `shapeConstr` would only be called within [cirls.fit][cirls.fit()], `intercept` is automatically determined from the [glm][stats::glm()] formula.
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
#' * [factor()]: for categorical variables. Extract the [contrasts][stats::contrasts()] to define the constraint matrix. here the `intercept` argument has the same interpretation as in the default method, i.e. if set to `TRUE` it means the `glm` model doesn't include an intercept externally to the factor. Note that, in this case, a simple dummy coding is done in R.
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
#' @returns A constraint matrix to be passed to `Cmat` in [cirls.fit][cirls.fit()].
#'
#' @references
#' Zhou, S. & Wolfe, D. A., 2000, On derivative estimation in spline regression. *Statistica Sinica* **10**, **93–108**.
#'
#' @seealso [cirls.fit()] which typically calls `shapeConstr` internally.
#'
#' @examples
#' # example code
#'
#'
#' @export
shapeConstr <- function(x, shape, ...) UseMethod("shapeConstr")
