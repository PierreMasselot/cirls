################################################################################
#
#  Initialize default parameters for coneproj
#
################################################################################

coneproj.def <- function(...){
  dots <- list(...)
  default <- list(msg = FALSE)
  utils::modifyList(default, dots)
}
