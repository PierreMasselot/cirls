################################################################################
#
# Method to create a constraint matrix for shape-constrained splines
# Works with various type of splines found in R.
#
################################################################################

#' Shape constraints
#'
#' @description
#' Builds a constraint matrix and bound vectors for defining shape constraints on a set of coefficients. This is a generic function designed to be used in the [constr][buildCmat()] formula interface. It allows methods for a wide range of regression terms (see Details).
#'
#' @param x An object representing a design matrix of predictor variables, typically basis functions. See Details for supported objects.
#' @param shape A character vector indicating one or several shape-constraints. See Details for supported shapes.
#' @param range A numeric vector of length 2. When specified, only this range is constrained. The default is to not restrict the range.
#' @param intercept For the default method, a logical value indicating if the design matrix includes an intercept. In most cases, it will be automatically extracted from `x`, but this argument can be used to override it.
#' @param ... Additional parameters passed to or from other methods.
#'
#' @details
#' This function is used to specify shape constraints on terms in a `cirls` model. Shapes can be enforced through the relation between the coefficients of several variables, for instance dummy indicators to impose shapes on the levels of a factor, or splines for smooth shapes. Note that this function can also be used to easily specify non-negativity or non-positivity constraints. See below for the list of implemented shapes and the examples section for several use cases.
#'
#' ## Usage
#'
#' The recommended usage is to use this function through a call to `shape` on a term in the [constr][buildCmat()] formula interface. This method is then called internally to create the constraint matrix and bound vectors. However, `shapeConstr` can also be called directly on a matrix-like object to manually build or inspect the constraint matrix.
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
#' ## Constraining a sub-range
#'
#' The `range` argument restricts the constraint over the specified range, leaving the association outside of this range unconstrained. Note that the exact range being constrained depends on the `x` and its specifications. For instance, when `x` represents bases that depend on knots or breakpoints (e.g. `bs`, `ns` or `strata`), the function will constrain between the knots/breakpoints closest to the specified range. In that sense `range` specifies a *minimum* range for the constraint, but most of the time a larger range will be constrained.
#'
#' ## Available methods
#'
#' In addition to the default method, `shapeConstr` currently supports methods for several classes, creating an appropriate shape-constraint matrix depending on the object. The full list (also provided by `methods(shapeConstr)`) is:
#'
#' ### Categorical variables
#' * [factor()]: for `factor` objects. Extract the [contrasts][stats::contrasts()] to define the constraint matrix. Here the `intercept` argument has the same interpretation as in the default method, i.e. if set to `TRUE` it means the `glm` model does not include an intercept externally to the factor. Note that, in this case, a simple dummy coding is done in R. Here `range` specifies the levels of constraining, having values between `1` and `nlevel(x)`.
#' * [strata][dlnm::strata()]: Indicator variables defining strata from the [dlnm][dlnm::dlnm] package. Here the shape is applied to the coefficients of the strata, treating it like a categorical variable. Here `range` is specified on the original scale of the variable, treating the strata similarly to spline bases.
#'
#' ### B-spline bases
#'
#' * [bs][splines::bs()]: B-splines bases from the [splines][splines()] package.
#' * [ns][splines::ns()]: Natural splines bases from the [splines][splines()] package.
#' * [ps][dlnm::ps()]: Penalised splines (P-Splines) from the [dlnm][dlnm::dlnm] package.
#'
#' ### Distributed-lag linear and nonlinear models (DLNM)
#'
#' * [onebasis][dlnm::onebasis()]: General method for basis functions generated in the package.
#' * [lin][dlnm::lin()]: Mostly for compatibility with the [dlnm][dlnm::dlnm] package. Here `"pos"` and `"inc"` (`"neg"` and `"dec"`) have the same interpretation, a non-negative (non-positive) coefficient associated to the linear term.
#' * [crossbasis][shapeConstr.crossbasis()]: Allows shape-constrained DLNM. See its dedicated [help page][shapeConstr.crossbasis()].
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub`, respectively).
#'
#' @seealso [buildCmat][buildCmat()] detailing the `constr` interface.
#' @aliases shape
#'
#' @references
#' Zhou, S. & Wolfe, D. A., 2000. On derivative estimation in spline regression. *Statistica Sinica* **10**, **93–108**.
#'
#' @example inst/examples/ex_london_nonneg.R
#' @example inst/examples/ex_warming_factor.R
#' @example inst/examples/ex_warming_splines.R
#' @example inst/examples/ex_warming_range.R
#'
#' @order 1
#' @export
shapeConstr <- function(x, ...) UseMethod("shapeConstr")
