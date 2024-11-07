################################################################################
#
#  Initialize default parameters for OSQP
#
################################################################################

# Also see: https://groups.google.com/g/osqp/c/BzEqWQR2dYY

osqp.def <- function(...){
  dots <- list(...)
  default <- list(verbose = FALSE, polish = TRUE)
  utils::modifyList(default, dots)
}
