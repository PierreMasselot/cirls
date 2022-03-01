#' @export
cirls.fit <- function (x, y, weights = rep.int(1, nobs), start = NULL,
  etastart = NULL, mustart = NULL, offset = rep.int(0, nobs),
  family = stats::gaussian(), control = list(), intercept = TRUE,
  singular.ok = TRUE)
{
  # Prepare CIRLS parameters
  control <- do.call("cirls.control", control)
  # Store variable names
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  # initialize convergence to FALSE (i.e. not converged)
  conv <- FALSE
  # Dimensions
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  # Initialize weights
  if (is.null(weights)) weights <- rep.int(1, nobs)
  # Initialize offset
  if (is.null(offset)) offset <- rep.int(0, nobs)
  # Intitalize family objects and check their validity
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)){
    stop("'family' argument seems not to be a valid family object",
      call. = FALSE)
  }
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if (is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  if (is.null(mustart)) {
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  # Set the results if there is no variables
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model",
        call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu))
      stop("invalid fitted means in empty model",
        call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    # If variables not empty, estimate model
    coefold <- NULL
    # Intialize eta by user provided value
    eta <- if (!is.null(etastart)){
      etastart
    } else {
      # Or using user provided starting coefficients
      if (!is.null(start)){
        if (length(start) != nvars){
          stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
            nvars, paste(deparse(xnames), collapse = ", ")),
            domain = NA)
        } else {
          coefold <- start
          offset + as.vector(if (NCOL(x) == 1L)
            x * start
            else x %*% start)
        }
        # if not use defaul starting values
      } else {
        family$linkfun(mustart)
      }
    }
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))){
      stop("cannot find valid starting values: please specify some",
        call. = FALSE)
    }
    # Initalize deviance for stopping criterion
    devold <- sum(dev.resids(y, mu, weights))
    # Initialize convergence flags
    boundary <- conv <- FALSE
    #### Loop for IRLS
    for (iter in 1L:control$maxit) {
      # Excluse observations with null weight
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (anyNA(varmu)) stop("NAs in V(mu)")
      if (any(varmu == 0)) stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good]))) stop("NAs in d(mu)/d(eta)")
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("no observations informative at iteration %d",
          iter), domain = NA)
        break
      }
      # Compute pseudo data
      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      # Compute pseudo weights
      w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
      # Weigh data
      wx <- w * x
      wz <- w * z
      ######################################################
      # PART SPECIFIC TO CONSTRAINED GLM
      ######################################################
      # Compute QR decomposition of design matrix
      wxqr <- qr(wx)
      Rmat <- qr.R(wxqr)
      effects <- qr.qty(wxqr, wz)
      # Pivoting in Cmat
      Cmat <- control$Cmat[,wxqr$pivot[seq_len(wxqr$rank)]]
      bvec <- control$bvec
      toremove <- apply(Cmat == 0, 1, all)
      if (any(toremove)){
        Cmat <- Cmat[!toremove]
        bvec <- bvec[!toremove]
      }
      # Fit QP
      fit <- quadprog::solve.QP(
        Dmat = crossprod(Rmat[seq_len(wxqr$rank),seq_len(wxqr$rank)]),
        dvec = crossprod(effects[seq_len(wxqr$rank)],
          Rmat[seq_len(wxqr$rank),seq_len(wxqr$rank)]),
        Amat = t(Cmat), bvec = bvec)
      # Check results
      if (any(!is.finite(fit$solution))) {
        conv <- FALSE
        warning(gettextf("non-finite coefficients at iteration %d",
          iter), domain = NA)
        break
      }
      # Check rank
      if (nobs < wxqr$rank) {
        stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation",
          "X matrix has rank %d, but only %d observations"),
          wxqr$rank, nobs), domain = NA)
      }
      if (!singular.ok && wxqr$rank < nvars) stop("singular fit encountered")
      # Update objects
      start <- rep(0, nvars)
      start[wxqr$pivot[seq_len(wxqr$rank)]] <- fit$solution
      ######################################################
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      # If required display advancement
      if (control$trace){
        cat("Deviance = ", dev, " Iterations - ",
          iter, "\n", sep = "")
      }
      boundary <- FALSE
      # If deviance is not finite, divide the step by two until it is finite
      if (!is.finite(dev)) {
        if (is.null(coefold)){
          stop("no valid set of coefficients has been found: please supply starting values",
            call. = FALSE)
        }
        warning("step size truncated due to divergence",
          call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit){
            stop("inner loop 1; cannot correct step size",
              call. = FALSE)
          }
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace){
          cat("Step halved: new deviance = ", dev,
            "\n", sep = "")
        }
      }
      # Idem, if valid coefficients not found, halve step
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold)){
          stop("no valid set of coefficients has been found: please supply starting values",
            call. = FALSE)
        }
        warning("step size truncated: out of bounds",
          call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit){
            stop("inner loop 2; cannot correct step size",
              call. = FALSE)
          }
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace){
          cat("Step halved: new deviance = ", dev,
            "\n", sep = "")
        }
      }
      # Check convergence
      if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      }
      else {
        devold <- dev
        coef <- coefold <- start
      }
    }
    # Give warnings for special cases
    if (!conv){
      warning("cirls.fit: algorithm did not converge",
        call. = FALSE)
    }
    if (boundary) {
      warning("cirls.fit: algorithm stopped at boundary value",
        call. = FALSE)
    }
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps)){
        warning("cirls.fit: fitted probabilities numerically 0 or 1 occurred",
          call. = FALSE)
      }
    }
    if (family$family == "poisson") {
      if (any(mu < eps)){
        warning("cirls.fit: fitted rates numerically 0 occurred",
          call. = FALSE)
      }
    }
    #######################################################################
    # Add warning if QR pivoting affects constraints (at last iteration)
    consrem <- apply(control$Cmat[,-wxqr$pivot[seq_len(wxqr$rank)]] != 0,
      2, any)
    if (any(consrem)){
      warning("some constraints removed because of rank deficiency",
        call. = FALSE)
    }
    # If final X not of full rank, assign NA coefficients
    if (wxqr$rank < nvars) coef[wxqr$pivot][seq.int(wxqr$rank + 1, nvars)] <- NA
    xxnames <- xnames[wxqr$pivot]
    #######################################################################
    # Residuals
    residuals <- (y - mu)/mu.eta(eta)
    # Rank
    wxqr$qr <- as.matrix(wxqr$qr)
    nr <- min(sum(good), nvars)
    # Update names
    names(coef) <- xnames
    colnames(wxqr$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  # Name to objects related to y
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  # Final pseudo weights including discarded observations
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  # names to weights
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY){
    names(effects) <- c(xxnames[seq_len(wxqr$rank)], rep.int("",
      sum(good) - wxqr$rank))
  }
  # Null deviance
  wtdmu <- if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  # Degrees of freedom
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) 0 else wxqr$rank
  resdf <- n.ok - rank
  aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
    effects = if (!EMPTY) effects, R = if (!EMPTY) Rmat, qr = if (!EMPTY) wxqr,
    rank = rank, family = family, linear.predictors = eta,
    deviance = dev, aic = aic.model, null.deviance = nulldev,
    iter = iter, inner.iter = fit$iterations[1],
    weights = wt, prior.weights = weights, df.residual = resdf,
    df.null = nulldf, y = y, converged = conv, boundary = boundary,
    ########################################################################
    coef.unconstrained = fit$unconstrained.solution, active.cons = fit$iact,
    Cmat = control$Cmat, bvec = control$bvec)
  ##########################################################################
}
