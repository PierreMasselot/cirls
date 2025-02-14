################################################################################
#
# Shape constraint matrix method:
# Method for onebasis
#
################################################################################

#' @rdname shapeConstr
#' @export
shapeConstr.onebasis <- function(x, shape, ...){

  # Extract the right method
  fun <- attr(x, "fun")
  met <- paste0("shapeConstr.", fun)
  if (!met %in% utils::methods("shapeConstr")) {
    warning(paste0("No existing 'shapeConstr' method for '", fun,
      "' functions. Using default method."))
    met <- "shapeConstr.default"
  }

  # Call the right method
  pars <- list(x = x, shape = shape, intercept = attr(x, "intercept"))
  pars <- utils::modifyList(pars, list(...))
  pars <- pars[names(pars) %in% names(formals(met))]
  do.call(met, pars)
}
