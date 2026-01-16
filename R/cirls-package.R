#' @description
#' The package `cirls` provides routines to fit **generalised linear models (GLM) with linear constraints on the coefficients**. These routines implement a general framework to perform a range of constraint-based methods, including common applications such as shape-constrained regression or compositional regression. However, the package allows users to specify *any* custom linear constraint for more niche applications. Functions in the `cirls` package include fitting procedures as well as dedicated methods to perform inference and model selection.
#'
#' The statistical framework is implemented by embedding constrained optimisation procedures in [glm][stats::glm()], so that users can flexibly perform constraint-based analyses with simple calls to this standard R function. The user can use the standard `glm` ecosystem on constrained model, augmented with dedicated methods for `cirls` objects.
#'
#' @section Methodological framework:
#' Applying linear constraints to GLMs entails restricting the fitted coefficients \eqn{\beta} to respect specific linear combinations. These constraints can represent assumptions on the coefficients or enforce identifiability. More formally, the fitted GLM is subject to:
#' \deqn{ \mathbf{l} \leq \mathbf{C\beta} \leq \mathbf{u} }
#' where \eqn{\mathbf{C}} is a matrix encoding the constraints, with each row representing a single constraint and each column representing a regressor in the design matrix of the model. \eqn{\mathbf{l}} and \eqn{\mathbf{u}} are lower and upper bounds for the constraints, respectively.
#'
#' The model is fitted by [constrained iteratively-reweighted least-squares][cirls.fit()] (CIRLS), which extends the classical IRLS algorithm to ensure the constraints specified by \eqn{\mathbf{C}}, \eqn{\mathbf{l}}, and \eqn{\mathbf{u}} are respected.
#'
#' Linearly constrained coefficients follow a **truncated multivariate Normal** (TMVN) distribution, *i.e.* a multivariate normal that respects the constraints specified by \eqn{\mathbf{C}}, \eqn{\mathbf{l}} and \eqn{\mathbf{u}}. Inference is conducted by [simulating][simulCoef()] from the TMVN distribution of fitted coefficients, and the output can then be used to compute useful summaries like the [variance-covariance matrix][vcov.cirls()] and [confidence intervals][confint.cirls()].
#'
#' @section Usage:
#'
#' ## Fitting the model
#' The `cirls` package has been developed as a **plug-in** to [glm][stats::glm()]. The main function is [cirls.fit][cirls.fit()], designed to be used within [glm][stats::glm()] through the argument `method = "cirls.fit"`. This will fit the GLM using specific CIRLS algorithms instead of the default IRLS algorithm. All arguments of [glm][stats::glm()] such as `formula`, `data` and `family` retain the same usage. Additional arguments specific to `cirls` are passed through the ellipsis `...` or the `control` argument (*but not both at the same time*). If any parameter is passed through `control`, then `...` will be ignored. See [cirls.control][cirls.control()] for the full list of additional arguments and examples below for how to use additional functions in the package.
#'
#' ## Specifying constraints
#' Constraints can be specified through arguments `constr` and/or `Cmat`/`lb`/`ub`. The [constr][buildCmat()] argument provides a user-friendly way to pass common built-in constraints through a formula and is the recommended way for new users. See the list of available built-in constraints below. Alternatively, the arguments [Cmat/lb/ub][buildCmat()] allow more flexibility in passing constraints. These arguments take either a full constraint matrix (`Cmat`) and bound vectors (`lb` and `ub`), or named lists to specify a constrained matrix and bound vectors for specific terms in the model. See the help for [buildCmat][buildCmat()] for more specific details.
#'
#' @section Functions included in the package:
#' [cirls.fit][cirls.fit()] is the workhorse performing the fitting, and is controlled by [cirls.control][cirls.control()]. The function returns an object of class `cirls` which inherits from class `glm`. This means that all methods available for usual `glm` objects can also be used on a `cirls` object. This typically includes [coef][stats::coef()] and [summary][stats::summary.glm()], for instance. However, some methods are supplanted by new methods that are specific to `cirls` objects. The sections below provides a comprehensive illustration of such specific methods, including an index of the functions.
#'
#' ## Inference
#' Inference for `cirls` objects includes computing the [variance-covariance matrix][vcov.cirls()] and [confidence intervals][confint.cirls()] for constrained (and unconstrained) coefficients. These functions build on [simulCoef][simulCoef()], which generates coefficient vectors based on the TMVN distribution of fitted coefficients, and can be used to compute other summaries.
#'
#' ## Model selection
#' Degrees of freedom of a fitted `cirls` object can be computed through the [edf][edf()] function. A specific [logLik][logLik.cirls()] method extracts the log-likelihood of the model with the right degrees of freedom in its `df` attribute. This can be used by [AIC][stats::AIC()] and [BIC][stats::BIC()] with appropriate degrees of freedom.
#'
#' ## Other functions
#' [uncons][uncons()] refits the model removing some or all of the constraints, acting similarly to the [update][stats::update()] function although focusing on the constraints.
#'
#' [buildCmat][buildCmat()] and [checkCmat][checkCmat()] are convenience functions that are mainly used internally in [cirls.fit][cirls.fit()]. They are not meant to be called directly, but can be useful to advanced users to build and/or check constraint matrices. Specifically, [buildCmat][buildCmat()] builds a full constraint matrix from the `constr` and/or `Cmat`/`lb`/`ub` arguments and a model frame. [checkCmat][checkCmat()] checks whether a constraint matrix can be reduced by detecting redundant and underlying equality constraints. The help pages of both functions also provide technical details on constraints and constraint matrices.
#'
#' @section Index:
#'
#' ## Built-in constraints
#' These constraints can be used in the `constr` formula as functions applied to specific terms:
#' * [shape][shapeConstr()]: constraining the shape of an association. Typically used with spline bases or factors.
#' * [zerosum][zerosumConstr()]: specify that a group of coefficients should sum to zero. Typically used for compositional regression.
#'
#' ## Methods specific to `cirls` objects
#'
#' The following supplant `glm` methods:
#' * [vcov][vcov.cirls()]: compute the variance-covariance matrix of coefficients.
#' * [confint][confint.cirls()]: compute confidence intervals for constrained coefficients.
#' * [logLik][logLik.cirls()]: extract the log-likelihood of a fitted object, assigning the correct [degrees of freedom][edf()].
#'
#' @section Datasets:
#'
#' * [fgl][fgl()]: Measurement of forensic glass fragments
#' * [london][london()]: Daily mortality, temperature, and pollution data in London
#' * [warming][warming()]: Global temperature anomaly
#'
#' @references
#' Masselot, P., Nenon, D., Vanoli, J., Chalabi, Z., Gasparrini, A., 2025. Estimation and inference in generalised linear models with constrained iteratively-reweighted least squares. *ArXiv preprint*. [DOI:10.48550/arXiv.2509.18406](https://doi.org/10.48550/arXiv.2509.18406)
#'
#' @example inst/examples/ex_warming_factor.R
#' @example inst/examples/ex_warming_splines.R
#' @example inst/examples/ex_fgl_coda.R
#' @example inst/examples/ex_london_nonneg.R
#'
"_PACKAGE"
