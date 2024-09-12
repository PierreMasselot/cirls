#' Constrained Iteratively Reweighted Least-Squares
#'
#' @description Fits a generalized linear model with linear constraints on the coefficients through a Constrained Iteratively Reweighted Least-Squares (CIRLS) algorithm.
#' This function is the constrained counterpart to [glm.fit][stats::glm.fit()] and is meant to be called by [glm][stats::glm()] through its `method` argument. See details for the main differences.
#'
#' @param x,y `x` is a design matrix and `y` is a vector of response observations. Usually internally computed by [glm][stats::glm()].
#' @param weights An optional vector of observation weights.
#' @param start Starting values for the parameters in the linear predictor.
#' @param etastart Starting values for the linear predictor.
#' @param mustart Starting values for the vector or means.
#' @param offset An optional vector specifying a known component in the model. See [model.offset][stats::model.offset()].
#' @param family The result of a call to a family function, describing the error distribution and link function of the model. See [family][stats::family()] for details of available family functions.
#' @param control A list of parameters controlling the fitting process. See details and [cirls.control][cirls.control()].
#' @param intercept Logical. Should an intercept be included in the null model?
#' @param singular.ok Logical. If `FALSE`, the function returns an error for singular fits.
#'
#' @details This function is a plug-in for [glm][stats::glm()] and works similarly to [glm.fit][stats::glm.fit()]. In addition to the parameters already available in [glm.fit][stats::glm.fit()], `cirls.fit` allows the specification of a constraint matrix `Cmat` with bound vectors `lb` and `ub` on the regression coefficients. These additional parameters can be passed through the `control` list or through `...` in [glm][stats::glm()].
#'
#' The CIRLS algorithm is a modification of the classical IRLS algorithm in which each update of the regression coefficients is performed by a quadratic program (QP), ensuring the update stays within the feasible region defined by `Cmat`, `lb` and `ub`. More specifically, this feasible region is defined as
#'
#' `lb <= Cmat %*% coefficients <= ub`
#'
#' where `coefficients` is the coefficient vector returned by the model. This specification allows for any linear constraint, including equality ones.
#'
#' ## Specifying `Cmat`, `lb` and `ub`
#'
#' `Cmat` is a matrix that defines the linear constraints. If provided directly as a matrix, the number of columns in `Cmat` must match the number of coefficients estimated by [glm][stats::glm()]. This includes all variables that are not involved in any constraint potential expansion such as factors or splines for instance, as well as the intercept. Columns not involved in any constraint will be filled by 0s.
#'
#' Alternatively, it may be more convenient to pass `Cmat` as a list of constraint matrices for specific terms. This is advantageous if a single term should be constrained in a model containing many terms. If provided as a list, `Cmat` is internally expanded to create the full constraint matrix. See examples of constraint matrices below.
#'
#' `lb` and `ub` are vectors defining the bounds of the constraints. By default they are set to `0` and `Inf`, meaning that the linear combinations defined by `Cmat` should be positive, but any bounds are possible. When some elements of `lb` and `ub` are identical, they define equality constraints. Setting `lb = -Inf` and `ub = Inf` disable the constraints.
#'
#' ## Quadratic programming solvers
#'
#' The function [cirls.fit][cirls.fit()] relies on a quadratic programming solver. Several solver are currently available.
#' - `"osqp"` (the default) solves the quadratic program via the Alternating Direction Method of Multipliers (ADMM). Internally it calls the function [solve_osqp][osqp::solve_osqp()].
#' - `"quadprog"` performs a dual algorithm to solve the quadratic program. It relies on the function [solve.QP][quadprog::solve.QP()].
#' - `"coneproj"` solves the quadratic program by a cone projection method. It relies on the function [qprog][coneproj::qprog()].
#'
#' Each solver has specific parameters that can be controlled through the argument `qp_pars`. Sensible defaults are set within [cirls.control][cirls.control()] and the user typically doesn't need to provide custom parameters.
#'
#' @return A `cirls` object inheriting from the class `glm`. At the moment, two non-standard methods specific to `cirls` objects are available: [vcov.cirls][vcov.cirls()] to obtain the coefficients variance-covariance matrix and [confint.cirls][confint.cirls()] to obtain confidence intervals. These custom methods account for the reduced degrees of freedom resulting from the constraints, see [vcov.cirls][vcov.cirls()] and [confint.cirls][confint.cirls()]. Any method for `glm` objects can be used, including the generic [coef][stats::coef()] or [summary][base::summary()] for instance.
#'
#' An object of class `cirls` includes all components from [glm][stats::glm()] objects, with the addition of:
#' \item{active.cons}{vector of indices of the active constraints in the fitted model.}
#' \item{inner.iter}{number of iterations performed by the last call to the QP solver.}
#' \item{Cmat, lb, ub}{the (expanded) constraint matrix, and lower and upper bound vectors.}
#'
#' @seealso [vcov.cirls][vcov.cirls()], [confint.cirls][confint.cirls()] for methods specific to `cirls` objects. [cirls.control][cirls.control()] for fitting parameters specific to [cirls.fit][cirls.fit()]. [glm][stats::glm()] for details on `glm` objects.
#'
#' @references
#' Goldfarb, D., Idnani, A., 1983. A numerically stable dual method for solving strictly convex quadratic programs. *Mathematical Programming* **27**, 1–33. \doi{10.1007/BF02591962}
#'
#' Meyer, M.C., 2013. A Simple New Algorithm for Quadratic Programming with Applications in Statistics. *Communications in Statistics - Simulation and Computation* **42**, 1126–1139. \doi{10.1080/03610918.2012.659820}
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd, S., 2020. OSQP: an operator splitting solver for quadratic programs. *Math. Prog. Comp.* **12**, 637–672. \doi{10.1007/s12532-020-00179-2}
#'
#' @example man/examples/cirls.fit.R
#'
#' @export
cirls.fit <- function (x, y, weights = rep.int(1, nobs), start = NULL,
  etastart = NULL, mustart = NULL, offset = rep.int(0, nobs),
  family = stats::gaussian(), control = list(), intercept = TRUE,
  singular.ok = TRUE)
{
  ################################################
  # Prepare CIRLS parameters
  control <- do.call("cirls.control", control)
  # Extract solver
  solver_fun <- sprintf("%s.fit", control$qp_solver)
  ################################################
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
      Cmat <- control$Cmat[,wxqr$pivot[seq_len(wxqr$rank)], drop = F]
      lb <- control$lb
      ub <- control$ub
      toremove <- apply(Cmat == 0, 1, all)
      if (any(toremove)){
        Cmat <- Cmat[!toremove,,drop = F]
        lb <- lb[!toremove]
        ub <- ub[!toremove]
      }
      # Fit QP
      fit <- do.call(solver_fun, list(
        Dmat = crossprod(Rmat[seq_len(wxqr$rank),seq_len(wxqr$rank)]),
        dvec = crossprod(effects[seq_len(wxqr$rank)],
          Rmat[seq_len(wxqr$rank),seq_len(wxqr$rank)]),
        Cmat = Cmat, lb = lb, ub = ub, qp_pars = control$qp_pars))
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
    consrem <- apply(control$Cmat[,-wxqr$pivot[seq_len(wxqr$rank)], drop = F] != 0,
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
  resdf <- n.ok - rank + length(fit$iact) # We remove the number of active constraints from the df of the model
  aic.model <- aic(y, nobs, mu, weights, dev) + 2 * rank
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
    effects = if (!EMPTY) effects, R = if (!EMPTY) Rmat, qr = if (!EMPTY) wxqr,
    rank = rank, family = family, linear.predictors = eta,
    deviance = dev, aic = aic.model, null.deviance = nulldev,
    iter = iter, weights = wt, prior.weights = weights, df.residual = resdf,
    df.null = nulldf, y = y, converged = conv, boundary = boundary,
    ########################################################################
    active.cons = fit$iact, inner.iter = fit$iterations,
    Cmat = control$Cmat, lb = control$lb, ub = control$ub, class = "cirls")
  ##########################################################################
}
