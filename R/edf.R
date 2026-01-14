################################################################################
#
#    Effective degrees of freedom
#
################################################################################

#' Expected degrees of freedom
#'
#' @description
#' Estimate expected degrees of freedom of a [cirls][cirls.fit()] object through simulations.
#'
#' @param object A `cirls` object or any object inheriting from `lm`, see details.
#' @param nsim The number of simulations.
#' @param seed An optional seed for the random number generator. See [set.seed][set.seed()].
#'
#' @details
#' Computes **unconstrained**, **observed** and **expected** degrees of freedom for `cirls` objects. The function also works for most objects inheriting from `lm` although in this case only the unconstrained (`udf`) makes sense.
#'
#' ## Unconstrained df
#'
#' **Unconstrained** degrees of freedom (`udf`) refer to the usual degrees of freedom, i.e. the number of estimated parameters \eqn{p}. In GLMs, this corresponds to the number of coefficients, plus one degree of freedom for the dispersion parameter when relevant.
#'
#' ## Observed df
#'
#' **Observed** degrees of freedom (`odf`) correspond to the constrained df of a fitted `cirls` model. This corresponds to \eqn{p-m_a}, where \eqn{m_a} is *the number of active constraint* in the fitted `cirls` model. Intuitively, the more restricted the final fit, the more degrees of freedom are decreased.
#'
#' ## Expected df
#'
#' **Expected** degrees of freedom (`edf`) account for the sampling variation in the number of active constraints. A different sample might result in a fit with another number of active constraints \eqn{m_a} which means `odf` represents an inaccurate estimation of the reduction in degrees of freedom induced by the constraints. `edf` is defined as \eqn{p - \bar{m_a}} where \eqn{\bar{m_a}} is the expected number of active constraints. Here, \eqn{\bar{m_a}} is estimated by sampling from the [unconstrained][uncons()] distribution of coefficients (see [simulCoef][simulCoef()]).
#'
#' @returns A vector of length three containing the three types of degrees of freedom:
#' \item{udf}{The *unconstrained* degrees of freedom, i.e. the rank plus any dispersion parameter for `glm` objects.}
#' \item{odf}{The *observed* degrees of freedom, that is `udf` minus the number of active constraints.}
#' \item{edf}{The *expected* degrees of freedom estimated by simulation as described in the details section.}
#' For `cirls` objects, the vector includes the simulated distribution of the number of active constraints as an `actfreq` attribute.
#'
#' @seealso [logLik.cirls][logLik.cirls()] which internally calls `edf` to compute degrees of freedom. [simulCoef][simulCoef()] for coefficient simulation.
#'
#' @references
#'  Meyer, M.C., 2013. Semi-parametric additive constrained regression. *Journal of Nonparametric Statistics* **25**, **715–730**. [DOI:10.1080/10485252.2013.797577](https://doi.org/10.1080/10485252.2013.797577)
#'
#' @example inst/examples/ex_warming_factor.R
#'
#' @export
edf <- function(object, nsim = 10000, seed = NULL){

  # Check object
  if (!inherits(object, "lm")) stop("'object' should inherit from 'lm'")

  # Extract unconstrained df (stats:::logLik.glm)
  fam <- stats::family(object)$family
  dispersion <- stats::family(object)$dispersion
  p <- object$rank
  if (!is.null(dispersion)) {
    if (is.na(dispersion))
      p <- p + 1
  } else if (fam %in% c("gaussian", "Gamma", "inverse.gaussian"))
    p <- p + 1

  # Start the result vector
  dfvec <- c(udf = p, odf = p - length(object$active.cons), edf = NA)

  #----- For cirls, now estimate reduced rank
  if (!is.null(object$Cmat)){

    # Simulate from the unconstrained model
    simures <- simulCoef(object, nsim = nsim, seed = seed, complete = TRUE,
      constrained = FALSE)

    # Check aliased coefficients
    aliased <- apply(is.na(simures), 2, all)
    if (all(aliased)){
      dfvec["edf"] <- NA
      return(dfvec)
    }

    # Remove aliased coefficients
    Cmat <- object$Cmat
    lb <- object$lb
    ub <- object$ub
    Cmat <- Cmat[,!aliased, drop = F]
    keep <- rowSums(Cmat != 0)
    Cmat <- Cmat[as.logical(keep),, drop = F]
    lb <- lb[as.logical(keep)]
    ub <- ub[as.logical(keep)]

    # Check active constraints
    cons <- Cmat %*% t(simures[,!aliased])
    active <- cons <= lb | cons >= ub

    # Compute the number of active constraints for each simulations
    actdist <- colSums(active)
    eact <- mean(actdist)

    # Compute average and store distribution
    dfvec["edf"] <- p - eact
    attr(dfvec, "actfreq") <- c(table(actdist)) / nsim
  }

  # Return
  dfvec
}
