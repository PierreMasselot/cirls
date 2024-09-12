################################################################################
#
#    Simulate coefficients from truncated multivariate normal
#
################################################################################

#' @rdname confint.cirls
#' @export
coef_simu <- function(object, nsim = 1000)
{
  # Extract original covariance matrix
  ovcov <- stats::summary.glm(object)$cov.scaled

  # Extract contraints
  Cmat <- object$Cmat
  lb <- object$lb
  ub <- object$ub

  # Check if any aliased coefficients and remove constraint if necessary
  aliased <- stats::summary.glm(object)$aliased
  betas <- stats::coef(object)[!aliased]
  Cmat <- Cmat[,!aliased, drop = F]
  keep <- rowSums(Cmat != 0)
  Cmat <- Cmat[as.logical(keep),, drop = F]
  lb <- lb[as.logical(keep)]
  ub <- ub[as.logical(keep)]

  #----- We simulate from a truncated normal with covariance CVCt

  # To allow back transformation, we "augment" the covariance matrix with
  #   its row null space
  Hmat <- nullspace(t(Cmat))
  Bmat <- rbind(Cmat, t(Hmat))

  # Transofrm V
  rectvcov <- Bmat %*% ovcov %*% t(Bmat)

  # Compute bounds of TMVN
  lbaug <- c(lb, rep(-Inf, ncol(Hmat)))
  lowervec <- lbaug - (Bmat %*% betas)
  ubaug <- c(ub, rep(Inf, ncol(Hmat)))
  uppervec <- ubaug - (Bmat %*% betas)

  # Initiate matrix with zeros for equality constraints
  truncres <- matrix(0, nrow = nsim, ncol = nrow(Bmat))
  eqind <- lowervec == uppervec

  # Simulate from truncated MVN
  truncres[,!eqind] <- suppressWarnings(TruncatedNormal::rtmvnorm(n = nsim,
    mu = rep(0, nrow(Bmat)), sigma = rectvcov, lb = lowervec, ub = uppervec))

  # Backtransform simulations
  backtruncres <- betas + solve(Bmat) %*% t(truncres)
  rownames(backtruncres) <- colnames(ovcov)

  # Export
  t(backtruncres)
}
