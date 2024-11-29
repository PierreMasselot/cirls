################################################################################
#
# Shape constraint matrix method:
# Method for onebasis
#
################################################################################

#' @export
cons_shape.onebasis <- function(x, diff = 0, sign = 1, ...){

  # Extract the right method
  fun <- attr(x, "fun")
  met <- paste0("cons_shape.", fun)
  if (!met %in% utils::methods("cons_shape")) {
    warning(paste0("No existing 'cons_shape' method for '", fun,
      "' functions. Using default method."))
    met <- "cons_shape.default"
  }

  # Call the right method
  pars <- list(x = x, diff = diff, sign = sign, intercept = attr(x, "intercept"))
  pars <- utils::modifyList(pars, list(...))
  pars <- pars[names(pars) %in% names(formals(met))]
  do.call(met, pars)
}
