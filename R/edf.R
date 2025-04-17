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
#' Simulates coefficient vectors from their unconstrained distribution, which is the non-truncated multivariate normal distribution. For each simulated vector, counts the number of violated constraints as the number of active constraints under the constrained distribution. The expected degrees of freedom is then the number of parameters minus the average number of active constraints.
#'
#' This procedure allows to account for the randomness of degrees of freedom for the constrained model. Indeed, the observed degrees of freedom is the number of parameters minus the number of active constraints. However, the number of active constraints is random as, some constraints can be active or not depending on the observed data. For instance, in a model for which the constraints are binding, the expected degrees of freedom will be close to the observed one, while in a model in which the constraints are irrelevant, the expected degrees of freedom will be closer to the unconstrained (usual) ones.
#'
#' # Note
#'
#' This function is implemented mainly for [cirls][cirls.fit()] objects and can return idiosyncratic results for other objects inheriting from `lm`. In this case, it will attempt to retrieve an 'edf' value, but simply return the rank of the model if this fails. For `glm` models for instance, it will return thrice the same value.
#'
#' @returns A vector of length three with components:
#' \item{udf}{The *unconstrained* degrees of freedom, i.e. the rank plus any dispersion parameter for `glm` objects.}
#' \item{odf}{The *observed* degrees of freedom, that is `udf` minus the number of active constraints.}
#' \item{edf}{The *expected* degrees of freedom estimated by simulation as described in the details section. For any other object inheriting from `lm`, attempts to retrieve the *effective* degrees of freedom.}
#' For `cirls` objects, the vector includes the simulated distribution of the number of active constraints as an `actfreq` attribute.
#'
#' @seealso [logLik.cirls][logLik.cirls()] which internally calls `edf` to compute degrees of freedom.
#'
#' @references
#'  Meyer, M.C., 2013. Semi-parametric additive constrained regression. *Journal of Nonparametric Statistics* **25**, **715–730**. \doi{10.1080/10485252.2013.797577}
#'
#' @example man/examples/edf_ex.R
#'
#' @export
edf <- function(object, nsim = 1000, seed = NULL){

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
  } else {
    #----- If not constrained, try and extract a value
    dfvec["edf"] <- object$edf %||% p
  }

  # Return
  dfvec
}
