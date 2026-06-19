################################################################################
# Function to apply constraints to crossbasis objects for DLNMs

#' Constraints for distributed lag linear and non-linear models
#'
#' @description
#' Methods to generate constraint matrices associated with [crossbasis][dlnm::crossbasis()] objects, allowing the fitting of constrained distributed lag linear and non-linear models (DLMs and DLNMs). Designed for use within the [constr][buildCmat()] interface in `cirls`.
#'
#' @param x A `crossbasis` object.
#' @param dim Either `"var"` (the default) or `"lag"`. This is the dimension on which the constraint will be applied, constraining either the exposure-response or lag-response relationship, respectively.
#' @param overall Only when `dim = "var"`, logical indicating whether the constraint should be applied only on the overall cumulative association or for each specific lag.
#' @param odrng A numeric vector of length 2 restricting the constraint on a specific range of the *other* dimension, namely the lag (for `dim = "var"`) or exposure (`dim = "lag"`) space. By default, no restriction is performed. See Details.
#' @param ... Parameters specific to the type of constraint to be applied. Includes for instance `shape` for [shapeConstr][shapeConstr()], or `value` for [boundConstr][boundConstr()]. See the main help page of the relevant generic method.
#'
#' @details
#' Constraints for DLMs and DLNMs apply to one of two dimensions of an exposure-lag-response association, specifically the exposure space (for `dim = "var"`) or the lag space (for `dim = "lag"`). Operationally, the functions generate a constraint matrix for the requested dimension of a `crossbasis` object, which is then expanded consistently in the other dimension. This allows the user to apply specific constraints to either the exposure-response association or the lag-response association, relying on the same [shape][shapeConstr()] or [bound][boundConstr()] methods used for other functions.
#'
#' ## Overall vs full surface
#'
#' By default, when `dim = "var"`, the constraint is applied across the whole lag space. This means that the exposure-response will be constrained at every lag. When `overall` is switched to `TRUE` (allowed only when `dim = "var"`), the constraint will only hold for the overall cumulative association (namely the net effect summed across lags), while it can be violated for exposure-responses defined at specific lags. See [crossreduce][dlnm::crossreduce()] for more information.
#'
#' ## Restricting the range on the other dimension
#'
#' The `odrng` argument can be used to restrict the constraint on a specific range of the other dimension. For instance, when `dim = "var"`, setting `odrng = c(0, 5)` will enforce the constraint on the exposure-response relationship only over lags 0 to 5, meaning it could be violated at other lags (> 5 in this example). Note that the actual range depends on the specific basis functions used for the opposite dimension (see the `range` argument in [shapeConstr][shapeConstr()] for additional details). When `overall = TRUE`, it will only constraint the cumulative association over the range specified by `odrng` (see above).
#'
#' ## Note
#'
#' These methods for [crossbasis][dlnm::crossbasis()] objects only allow to specify a single constraint to one dimension of the exposure-lag-response association. However, the user can apply multiple constraints, possibly to both dimensions, by adding multiple terms for the same object in the [constr][buildCmat()] formula (see Examples below).
#'
#' ## Warnings
#'
#' In the current implementation, constraints applied to the lag dimension (when `dim = "lag"`) can be straightforwardly defined when the exposure-response function is linear, namely for simpler DLMs. For non-linear exposure-responses, the centering performed by [crosspred][dlnm::crosspred()] when predicting the exposure-lag-response relationship can break the constraint definition and return unpredictable results. Users are advised not to apply constraints in the lag dimension of complex DLNMs, or to do so with particular care.
#'
#' @returns A list containing the constraint matrix `Cmat`, and lower/upper bound vectors (`lb` and `ub`, respectively).
#'
#' @seealso [crossbasis][dlnm::crossbasis()] for the definition of DLMs and DLNMs as well as their dimensions. [shapeConstr][shapeConstr()] and [boundConstr][boundConstr()] for generic constraint methods. [buildCmat][buildCmat()] for how to specify constraints.
#'
#' @references
#' Gasparrini, A., Armstrong, B., Kenward, M.G., 2010. Distributed lag non-linear models. *Statistics in Medicine* **29**, **2224–2234**. [DOI:10.1002/sim.3940](https://doi.org/10.1002/sim.3940)
#'
#' Gasparrini, A., Armstrong, B., 2013. Reducing and meta-analysing estimates from distributed lag non-linear models. *BMC Medical Research Methodology* **13**, **1**. [DOI:10.1186/1471-2288-13-1](https://doi.org/10.1186/1471-2288-13-1)
#'
#' @example inst/examples/ex_london_dlm.R
#' @example inst/examples/ex_london_dlnm.R
#'
#' @export
shapeConstr.crossbasis <- function(x, dim = "var", overall = FALSE, odrng = NULL,
  ...)
{
  # Call cbConstr
  cbConstr(x, constr = "shape", pars = list(...), dim = dim, overall = overall,
    odrng = odrng)
}
