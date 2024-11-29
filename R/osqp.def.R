################################################################################
#
#  Initialize default parameters for OSQP
#
################################################################################

# Also see: https://groups.google.com/g/osqp/c/BzEqWQR2dYY

osqp.def <- function(...){
  dots <- list(...)
  default <- list(verbose = FALSE, polish = TRUE,
    eps_abs = 1e-6, eps_rel = 1e-6)
  utils::modifyList(default, dots)
}
