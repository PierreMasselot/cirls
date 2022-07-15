\name{cirls.control}
\alias{cirls.control}

\title{Parameters controlling the cirls fit}
\description{
Internal function setting the parameters to control the cirls fit. It is called internally by \code{\link{cirls.fit}} and should not be called directly by the user. Can be called through the \code{control} argument or through \code{...} in \code{\link[stats]{glm}}.
}

\usage{
cirls.control <- function (epsilon = 1e-08, maxit = 25, trace = FALSE,
  Cmat = NULL, lb = 0L, ub = Inf, qp_solver = "osqp", qp_pars = list())
}
