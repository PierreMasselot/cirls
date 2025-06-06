% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cirls.fit.R
\name{cirls.fit}
\alias{cirls.fit}
\title{Constrained Iteratively Reweighted Least-Squares}
\usage{
cirls.fit(x, y, weights = rep.int(1, nobs), start = NULL,
  etastart = NULL, mustart = NULL, offset = rep.int(0, nobs),
  family = stats::gaussian(), control = list(), intercept = TRUE,
  singular.ok = TRUE)
}
\arguments{
\item{x, y}{\code{x} is a design matrix and \code{y} is a vector of response observations. Usually internally computed by \link[stats:glm]{glm}.}

\item{weights}{An optional vector of observation weights.}

\item{start}{Starting values for the parameters in the linear predictor.}

\item{etastart}{Starting values for the linear predictor.}

\item{mustart}{Starting values for the vector or means.}

\item{offset}{An optional vector specifying a known component in the model. See \link[stats:model.extract]{model.offset}.}

\item{family}{The result of a call to a family function, describing the error distribution and link function of the model. See \link[stats:family]{family} for details of available family functions.}

\item{control}{A list of parameters controlling the fitting process. See details and \link[=cirls.control]{cirls.control}.}

\item{intercept}{Logical. Should an intercept be included in the null model?}

\item{singular.ok}{Logical. If \code{FALSE}, the function returns an error for singular fits.}
}
\value{
A \code{cirls} object inheriting from the class \code{glm}. At the moment, two non-standard methods specific to \code{cirls} objects are available: \link[=vcov.cirls]{vcov.cirls} to obtain the coefficients variance-covariance matrix and \link[=confint.cirls]{confint.cirls} to obtain confidence intervals. These custom methods account for the reduced degrees of freedom resulting from the constraints, see \link[=vcov.cirls]{vcov.cirls} and \link[=confint.cirls]{confint.cirls}. Any method for \code{glm} objects can be used, including the generic \link[stats:coef]{coef} or \link[base:summary]{summary} for instance.

An object of class \code{cirls} includes all components from \link[stats:glm]{glm} objects, with the addition of:
\item{Cmat, lb, ub}{the constraint matrix, and lower and upper bound vectors. If provided as lists, the full expanded matrix and vectors are returned.}
\item{active.cons}{vector of indices of the active constraints in the fitted model.}
\item{inner.iter}{number of iterations performed by the last call to the QP solver.}
\item{etastart}{the initialisation of the linear predictor \code{eta}. The same as \code{etastart} when provided.}
\item{singular.ok}{the value of the \code{singular.ok} argument.}
}
\description{
Fits a generalized linear model with linear constraints on the coefficients through a Constrained Iteratively Reweighted Least-Squares (CIRLS) algorithm.
This function is the constrained counterpart to \link[stats:glm]{glm.fit} and is meant to be called by \link[stats:glm]{glm} through its \code{method} argument. See details for the main differences.
}
\details{
This function is a plug-in for \link[stats:glm]{glm} and works similarly to \link[stats:glm]{glm.fit}. In addition to the parameters already available in \link[stats:glm]{glm.fit}, \code{cirls.fit} allows the specification of a constraint matrix \code{Cmat} with bound vectors \code{lb} and \code{ub} on the regression coefficients. These additional parameters can be passed through the \code{control} list or through \code{...} in \link[stats:glm]{glm} \emph{but not both}. If any parameter is passed through \code{control}, then \code{...} will be ignored.

The CIRLS algorithm is a modification of the classical IRLS algorithm in which each update of the regression coefficients is performed by a quadratic program (QP), ensuring the update stays within the feasible region defined by \code{Cmat}, \code{lb} and \code{ub}. More specifically, this feasible region is defined as

\verb{lb <= Cmat \%*\% coefficients <= ub}

where \code{coefficients} is the coefficient vector returned by the model. This specification allows for any linear constraint, including equality ones.
\subsection{Specifying constraints}{

The package includes several mechanisms to specify constraints. The most straightforward is to pass a matrix to \code{Cmat} with associated bound vectors in \code{lb} and \code{ub}. In this case, the number of columns in \code{Cmat} must match the number of coefficients estimated by \link[stats:glm]{glm}. This includes all variables that are not involved in any constraint potential expansion such as factors or splines for instance, as well as the intercept. By default \code{lb} and \code{ub} are set to \code{0} and \code{Inf}, respectively, meaning that the linear combinations defined by \code{Cmat} should be positive, but any bounds are possible. When some elements of \code{lb} and \code{ub} are identical, they define equality constraints. Setting \code{lb = -Inf} and \code{ub = Inf} disable the constraints.

In many cases, it is however more convenient to use the argument \code{constr} which allow specifying constraints through a formula. Additionally, \code{Cmat} (as well as \code{lb} and \code{ub}) can be passed as lists of matrices for specific terms. See \link[=buildCmat]{buildCmat} for full details on how to specify constraints.
}

\subsection{Quadratic programming solvers}{

The function \link[=cirls.fit]{cirls.fit} relies on a quadratic programming solver. Several solver are currently available.
\itemize{
\item \code{"quadprog"} (the default) performs a dual algorithm to solve the quadratic program. It relies on the function \link[quadprog:solve.QP]{solve.QP}.
\item \code{"osqp"} solves the quadratic program via the Alternating Direction Method of Multipliers (ADMM). Internally it calls the function \link[osqp:solve_osqp]{solve_osqp}.
\item \code{"coneproj"} solves the quadratic program by a cone projection method. It relies on the function \link[coneproj:qprog]{qprog}.
}

Each solver has specific parameters that can be controlled through the argument \code{qp_pars}. Sensible defaults are set within \link[=cirls.control]{cirls.control} and the user typically doesn't need to provide custom parameters. \code{"quadprog"} is set as the default being generally more reliable than the other solvers. \code{"osqp"} is faster but can be less accurate, in which case it is recommended to increase convergence tolerance at the cost of speed.
}
}
\examples{
####################################################
# Simple non-negative least squares

# Simulate predictors and response with some negative coefficients
set.seed(111)
n <- 100
p <- 10
betas <- rep_len(c(1, -1), p)
x <- matrix(rnorm(n * p), nrow = n)
y <- x \%*\% betas + rnorm(n)

# Define constraint matrix (includes intercept)
# By default, bounds are 0 and +Inf
Cmat <- cbind(0, diag(p))

# Fit GLM by CIRLS
res1 <- glm(y ~ x, method = cirls.fit, Cmat = Cmat)
coef(res1)

# Same as passing Cmat through the control argument
res2 <- glm(y ~ x, method = cirls.fit, control = list(Cmat = Cmat))
identical(coef(res1), coef(res2))

####################################################
# Increasing coefficients

# Generate two group of variables: an isotonic one and an unconstrained one
set.seed(222)
p1 <- 5; p2 <- 3
x1 <- matrix(rnorm(100 * p1), 100, p1)
x2 <- matrix(rnorm(100 * p2), 100, p2)

# Generate coefficients: those in b1 should be increasing
b1 <- runif(p1) |> sort()
b2 <- runif(p2)

# Generate full data
y <- x1 \%*\% b1 + x2 \%*\% b2 + rnorm(100, sd = 2)

#----- Fit model

# Create constraint matrix and expand for intercept and unconstrained variables
Ciso <- diff(diag(p1))
Cmat <- cbind(0, Ciso, matrix(0, nrow(Ciso), p2))

# Fit model
resiso <- glm(y ~ x1 + x2, method = cirls.fit, Cmat = Cmat)
coef(resiso)

# Compare with unconstrained
plot(c(0, b1, b2), pch = 16)
points(coef(resiso), pch = 16, col = 3)
points(coef(glm(y ~ x1 + x2)), col = 2)

#----- More convenient specification

# Cmat can be provided as a list
resiso2 <- glm(y ~ x1 + x2, method = cirls.fit, Cmat = list(x1 = Ciso))

# Internally Cmat is expanded and we obtain the same result
identical(resiso$Cmat, resiso2$Cmat)
identical(coef(resiso), coef(resiso2))

#----- Adding bounds to the constraints
# Difference between coefficients must be above a lower bound and below 1
lb <- 1 / (p1 * 2)
ub <- 1

# Re-fit the model
resiso3 <- glm(y ~ x1 + x2, method = cirls.fit, Cmat = list(x1 = Ciso),
  lb = lb, ub = ub)

# Compare the fit
plot(c(0, b1, b2), pch = 16)
points(coef(resiso), pch = 16, col = 3)
points(coef(glm(y ~ x1 + x2)), col = 2)
points(coef(resiso3), pch = 16, col = 4)
}
\references{
Goldfarb, D., Idnani, A., 1983. A numerically stable dual method for solving strictly convex quadratic programs. \emph{Mathematical Programming} \strong{27}, 1–33. \doi{10.1007/BF02591962}

Meyer, M.C., 2013. A Simple New Algorithm for Quadratic Programming with Applications in Statistics. \emph{Communications in Statistics - Simulation and Computation} \strong{42}, 1126–1139. \doi{10.1080/03610918.2012.659820}

Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd, S., 2020. OSQP: an operator splitting solver for quadratic programs. \emph{Math. Prog. Comp.} \strong{12}, 637–672. \doi{10.1007/s12532-020-00179-2}
}
\seealso{
\link[=vcov.cirls]{vcov.cirls}, \link[=confint.cirls]{confint.cirls} for methods specific to \code{cirls} objects. \link[=cirls.control]{cirls.control} for fitting parameters specific to \link[=cirls.fit]{cirls.fit}. \link[stats:glm]{glm} for details on \code{glm} objects.
}
