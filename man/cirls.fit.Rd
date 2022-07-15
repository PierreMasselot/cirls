\name{cirls.fit}
\alias{cirls.fit}

\title{Constrained Iteratively Reweighted Least-Squares}
\description{
Fits a generalized linear model with linear constraints on the coefficients through a constrained iteratively reweighted least-squares algorithm. 
Plug-in function meant to be called by \code{\link[stats]{glm}} through its \code{method} argument.
}


\usage{
cirls.fit <- function (x, y, weights = rep.int(1, nobs), start = NULL,
  etastart = NULL, mustart = NULL, offset = rep.int(0, nobs),
  family = stats::gaussian(), control = list(), intercept = TRUE,
  singular.ok = TRUE)
}
