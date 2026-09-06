################################################################################
#
#    Simulate coefficients from truncated multivariate normal
#
################################################################################

#' Methods for inference on the coefficients of a `cirls` object.
#'
#' @description
#' `simulCoef` simulates coefficients for a fitted `cirls` object and uses these simulations for inference. `confint` and `vcov` directly compute confidence intervals and the variance-covariance matrix for coefficients from a fitted `cirls` object. These methods for `cirls` objects supersede the default `glm` methods.
#'
#' @param object An object of class `cirls` or `sim.cirls`.
#' @param nsim The number of simulations to perform.
#' @param seed Either `NULL` or an integer that will be used in a call to [set.seed()] before simulating the coefficients.
#' @param complete If `FALSE`, it does not return inference for undetermined coefficients in case of an over-determined model.
#' @param parm A specification of which parameters to compute the confidence intervals for. Either a vector of index numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required.
#' @param constrained A logical switch indicating whether to simulate from the constrained (the default) or unconstrained (see [uncons][uncons()]) coefficients distribution. When set to `FALSE` in `vcov`, returns the untruncated covariance matrix as in an unconstrained GLM.
#' @param ... Further arguments passed to or from other methods. For `vcov` and `confint`, it can be used to provide a `seed` for the internal coefficient simulation.
#'
#' @details
#'
#' ## Coefficient inference
#' To perform inference for coefficients the `simulCoef` function simulates from the distribution of \eqn{\mathbf{C}\beta} which follows a **Truncated Multivariate Normal** distribution \eqn{TMVN(\mathbf{C}\beta^{*}, \mathbf{C}\mathbf{\Sigma}^{*}\mathbf{C}^{T}, \mathbf{l}, \mathbf{u})} where \eqn{\mathbf{C}} is the constraint matrix with bound vectors \eqn{\mathbf{l}} and \eqn{\mathbf{u}}, and \eqn{\beta^{*}} and \eqn{\mathbf{\Sigma}^{*}} are the unconstrained coefficient vector and variance matrix. The TMVN simulations are then back-transformed to the domain of \eqn{\beta} to allow for inference.
#'
#' ## Functions
#' `simulCoef` is the workhorse of the inference, performing the simulations described above. It is called internally by `confint.cirls` and `vcov.cirls` to compute confidence intervals and the variance-covariance matrix, respectively. These are custom methods for [cirls][cirls.fit()] objects that supersede the default methods used for [glm][stats::glm()] objects. Alternatively `simulCoef` can be used directly to simulate `nsim` coefficient to then be passed to `vcov.sim.cirls` and `confint.sim.cirls` methods. This avoids simulating several times when both confidence intervals and variance-covariance are needed for instance.
#'
#' @note
#' These methods only work when there are less constraints than variables in `cirls` model, i.e. when `Cmat` has less rows than columns.
#'
#' @returns
#' For `simulCoef`: a `sim.cirls` object, that is a matrix of `nsim` rows containing simulated coefficients with the attributes `seed`, `complete` and `constrained`.
#'
#' For `confint` methods: a two-column matrix with columns giving lower and upper confidence limits for each parameter.
#'
#' For `vcov` methods: a matrix of the estimated covariance between the parameter estimates of the model.
#'
#' @references
#' Geweke, J.F., 1996. Bayesian Inference for Linear Models Subject to Linear Inequality Constraints, in: Lee, J.C., Johnson, W.O., Zellner, A. (Eds.), Modelling and Prediction Honoring Seymour Geisser. *Springer, New York, NY*, pp. 248–263. [DOI:10.1007/978-1-4612-2414-3_15](https://doi.org/10.1007/978-1-4612-2414-3_15)
#'
#' Botev, Z.I., 2017, The normal law under linear restrictions: simulation and estimation via minimax tilting, *Journal of the Royal Statistical Society, Series B*, **79** (**1**), pp. 1–24. [DOI:10.1111/rssb.12162](https://doi.org/10.1111/rssb.12162)
#'
#' @seealso [rtmvnorm][TruncatedNormal::tmvnorm()] which is the routine used internally to simulate from a TMVN. [reduceCons][reduceCons()] to reduce the set of constraints.
#'
#' @example inst/examples/ex_warming_inference.R
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

    # Sometimes, simvcov is not symmetric due to small numerical differences
    if(!isSymmetric(simvcov)) simvcov <- (simvcov + t(simvcov)) / 2

    # Expand bounds for simulation
    lowervec <- c(lb, rep(-Inf, ncol(Hmat)))
    uppervec <- c(ub, rep(Inf, ncol(Hmat)))

    # Initiate matrix that contains results
    truncres <- matrix(NA, nrow = nsim, ncol = nrow(Bmat))

    # Fill up the equality constrained variables
    eqind <- (uppervec - lowervec) < sqrt(.Machine$double.eps)
    truncres[,eqind] <- rep(lowervec[eqind], each = nsim)

    # If there are equality constraints, use Schur complement
    if (any(eqind)){
      simbeta <- c(simbeta[!eqind] + simvcov[!eqind, eqind, drop = FALSE] %*%
        solve(simvcov[eqind, eqind, drop = FALSE]) %*%
        (lb[eqind] - simbeta[eqind]))
      simvcov <- simvcov[!eqind, !eqind, drop = FALSE] -
        simvcov[!eqind, eqind, drop = FALSE] %*%
        solve(simvcov[eqind, eqind, drop = FALSE]) %*%
        simvcov[eqind, !eqind, drop = FALSE]
    }

    # Simulate from truncated MVN, use the internal function directly
    # If a warning is thrown, stop the simulation and return NAs
    tmvnsim <- tryCatch(TruncatedNormal::mvrandn(n = nsim,
        l = lowervec[!eqind], u = uppervec[!eqind],
        Sig = simvcov, mu = simbeta),
      warning = function(w){
        warning(paste0("Could not sample from TMVN, ",
          "possibly due to instability in the unconstrained model"))
        matrix(NA, length(simbeta), nsim)
      })
    truncres[, !eqind] <- t(tmvnsim)

    # truncres[,!eqind, drop = FALSE] <- suppressWarnings(
    #   TruncatedNormal::rtmvnorm(n = nsim, mu = simbeta, sigma = simvcov,
    #     lb = lowervec, ub = uppervec, check = FALSE)
    # )

    # Backtransform simulations
    simu <- t(solve(Bmat) %*% t(truncres))
  } else {

    #----- If not, simulate without bounds
    simu <- TruncatedNormal::rtmvnorm(n = nsim, mu = ubeta, sigma = uvcov,
      check = FALSE)
  }

  #----- Return, including NAs if complete == TRUE

  # Add NAs if aliased coefficients
  if (complete) {
    outsimu <- matrix(NA, nrow = nsim, ncol = length(aliased),
      dimnames = list(NULL, names(aliased)))
    outsimu[, which(!aliased)] <- simu
  } else {
    outsimu <- simu
    colnames(outsimu) <- names(aliased[!aliased])
  }

  # Object
  attributes(outsimu) <- c(attributes(outsimu), list(class = "sim.cirls",
    seed = seed, complete = complete, constrained = constrained))

  # Export
  outsimu
}
