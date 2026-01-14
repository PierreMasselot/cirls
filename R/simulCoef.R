################################################################################
#
#    Simulate coefficients from truncated multivariate normal
#
################################################################################

#' Methods for inference on the coefficients of a `cirls` object.
#'
#' @description
#' Simulates coefficients for a fitted `cirls` object and uses these simulations for inference. `confint` and `vcov` directly compute confidence intervals and the Variance-Covariance matrix for coefficients from a fitted `cirls` object. These methods supersede the default methods for `cirls` objects.
#'
#' @param object A fitted `cirls` object.
#' @param nsim The number of simulations to perform.
#' @param seed Either NULL or an integer that will be used in a call to [set.seed()] before simulating the coefficients.
#' @param complete If FALSE, doesn't return inference for undetermined coefficients in case of an over-determined model.
#' @param parm A specification of which parameters to compute the confidence intervals for. Either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required.
#' @param constrained A logical switch indicating whether to simulate from the constrained (the default) or unconstrained (see [uncons][uncons()]) coefficients distribution. When set to `FALSE` in `vcov`, returns the untruncated covariance matrix as in an unconstrained GLM.
#' @param ... Further arguments passed to or from other methods. For `vcov` and `confint` can be used to provide a `seed` for the internal coefficient simulation.
#'
#' @details
#'
#' ## Coefficient inference
#' To perform inference for coefficients the `simulCoef` function simulates from the distribution of \eqn{\mathbf{C}\beta} which follows a **Truncated Multivariate Normal** distribution \eqn{TMVN(\mathbf{C}\beta^{*}, \mathbf{C}\mathbf{\Sigma}^{*}\mathbf{C}^{T}, \mathbf{l}, \mathbf{u})} where \eqn{\mathbf{C}} is the constraint matrix with bound vectors \eqn{\mathbf{l}} and \eqn{\mathbf{u}}, and \eqn{\beta^{*}} and \eqn{\mathbf{\Sigma}^{*}} are the unconstrained coefficient vector and variance matrix. The TMVN simulations are then back-transformed to the domain of \eqn{\beta} to allow for inference.
#'
#' ## Functions
#' `simulCoef` is the workhorse of the inference and is called internally by `confint` and `vcov`. They are custom methods for [cirls][cirls.fit()] objects to supersede the default methods used for [glm][stats::glm()] objects. `simulCoef` doesn't need to be used directly for confidence intervals and variance-covariance matrices, but it can be used to check other summaries of the coefficients distribution.
#'
#' @note
#' These methods only work when there are less constraints than variables in `cirls` model, i.e. when `Cmat` has less rows than columns.
#'
#' @returns
#' For `simulCoef`, a matrix with `nsim` rows containing simulated coefficients.
#'
#' For `confint`, a two-column matrix with columns giving lower and upper confidence limits for each parameter.
#'
#' For `vcov`, a matrix of the estimated covariances between the parameter estimates of the model.
#'
#' @references
#' Geweke, J.F., 1996. Bayesian Inference for Linear Models Subject to Linear Inequality Constraints, in: Lee, J.C., Johnson, W.O., Zellner, A. (Eds.), Modelling and Prediction Honoring Seymour Geisser. *Springer, New York, NY*, pp. 248–263. [DOI:10.1007/978-1-4612-2414-3_15](https://doi.org/10.1007/978-1-4612-2414-3_15)
#'
#' Botev, Z.I., 2017, The normal law under linear restrictions: simulation and estimation via minimax tilting, *Journal of the Royal Statistical Society, Series B*, **79** (**1**), pp. 1–24. [DOI:10.1111/rssb.12162](https://doi.org/10.1111/rssb.12162)
#'
#' @seealso [rtmvnorm][TruncatedNormal::tmvnorm()] which is the routine used internally to simulate from a TMVN. [checkCmat][checkCmat()] to check if the constraint matrix can be reduced.
#'
#' @example inst/examples/ex_warming_factor.R
#'
#' @order 1
#' @export
simulCoef <- function(object, nsim = 1, seed = NULL, complete = TRUE,
  constrained = TRUE)
{

  # Extract unconstrained coefficients
  ufit <- uncons(object)
  ubeta <- stats::coef(ufit, complete = FALSE)
  aliased <- stats::summary.glm(ufit)$aliased
  uvcov <- stats::.vcov.aliased(aliased, stats::summary.glm(ufit)$cov.scaled,
    complete = FALSE)

  # Check uvcov exists
  if (any(is.na(uvcov))){
    warning("Impossible to perform inference: unconstrained vcov matrix undefined. Returning NAs")
    return(matrix(NA, nsim, ifelse(complete, length(aliased), sum(aliased)),
      dimnames = list(NULL, names(aliased))))
  }

  # Set seed if necessary
  if (!is.null(seed)){
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  # Extract contraints
  Cmat <- object$Cmat
  lb <- object$lb
  ub <- object$ub

  # If Cmat is empty, switch of the constrained simulation
  if (NROW(Cmat) == 0) constrained <- FALSE

  #----- Extract constraints and transform to "square" domain then simulate
  if (constrained){

    # Remove constraints affected by aliased coefficients
    Cmat <- Cmat[,!aliased, drop = F]
    keep <- rowSums(Cmat != 0)
    Cmat <- Cmat[as.logical(keep),, drop = F]
    lb <- lb[as.logical(keep)]
    ub <- ub[as.logical(keep)]

    # Check constraint matrix
    rowrk <- qr(t(Cmat))$rank
    if (nrow(Cmat) > rowrk){
      warning("Impossible to perform inference: constraint matrix not of full row rank. Returning NAs")
      return(matrix(NA, nsim, ifelse(complete, length(aliased), sum(aliased)),
        dimnames = list(NULL, names(aliased))))
    }

    # To allow back transformation, we "augment" the constraint matrix with
    #   its row null space (Tallis 1965)
    Hmat <- nullspace(t(Cmat))
    Bmat <- rbind(Cmat, t(Hmat))

    # Transform parameters
    simvcov <- Bmat %*% uvcov %*% t(Bmat)
    simbeta <- drop(Bmat %*% ubeta)

    # Expand bounds for simulation
    lowervec <- c(lb, rep(-Inf, ncol(Hmat)))
    uppervec <- c(ub, rep(Inf, ncol(Hmat)))

    # Initiate matrix with zeros for equality constraints
    truncres <- matrix(0, nrow = nsim, ncol = nrow(Bmat))
    eqind <- lowervec == uppervec

    # Simulate from truncated MVN
    truncres[,!eqind] <- suppressWarnings(TruncatedNormal::rtmvnorm(n = nsim,
      mu = simbeta, sigma = simvcov, lb = lowervec, ub = uppervec,
      check = FALSE))

    # Backtransform simulations
    simu <- t(solve(Bmat) %*% t(truncres))
  } else {

    #----- If not, simulate without bounds
    simu <- TruncatedNormal::rtmvnorm(n = nsim, mu = ubeta, sigma = uvcov,
      check = FALSE)
  }

  #----- Return, including NAs if complete == TRUE

  if (complete) {
    outsimu <- matrix(NA, nrow = nsim, ncol = length(aliased),
      dimnames = list(NULL, names(aliased)))
    outsimu[, which(!aliased)] <- simu
  } else {
    outsimu <- simu
    colnames(outsimu) <- names(aliased[!aliased])
  }

  # Export
  outsimu
}
