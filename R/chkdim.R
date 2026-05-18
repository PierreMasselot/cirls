################################################################################
#
# Check dimension for constraints on crossbases
#
################################################################################

# dim is the character providing the dimension
# dfs represent the dimensions of the crossbasis

chkdim <- function(dim = NULL, dfs){

  # By default, if var is linear, then apply constraint on lag
  # Removed as can cause problems with nonlinear bases with 1df
  # if (is.null(dim)){
  #   dim <- if(dfs[1] == 1) "lag" else "var"
  # }

  # Check it is admissible
  dim <- match.arg(dim, c("var", "lag"))


  # Warning message for lag constraint in DLNM
  if (dim == "lag" && dfs[1] > 1) warning(paste0("When var is nonlinear, ",
    "constraints on lag can be broken by centering"))

  # Return
  dim
}
