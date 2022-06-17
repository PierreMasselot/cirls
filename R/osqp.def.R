################################################################################
#
#  Initialize default parameters for OSQP
#
################################################################################

osqp.def <- function(...){
  dots <- list(...)
  default <- list(verbose = FALSE, polish = TRUE)
  utils::modifyList(default, dots)
}
