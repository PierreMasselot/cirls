################################################################################
#
# Function to create constraint matrices for cross-bases
# Used as internal function for other methods
#
################################################################################

# General function to apply a constraint on one dimension of a crossbasis
# x is the crossbasis
# constr is the constr function to call
# pars include the parameters for the constr function
# dim is the dimension on which to apply the constraint
# overall refers to whether the constraint has to be applied on the overall
#     cumulative or everywhere
# odrng to provide a specific range (of the other dimension) on which to apply
#     the constraint

cbConstr <- function(x, constr, pars = list(), dim = "var", overall = FALSE,
  odrng = NULL)
{

  # Get info from crossbasis
  cbattr <- attributes(x)
  varrng <- cbattr$range
  lagrng <- cbattr$lag
  dfs <- cbattr$df

  # Check dim parameter
  dim <- chkdim(dim, dfs)

  # Template basis for var (it won't change anyway)
  varbasis <- do.call(dlnm::onebasis, c(list(x = varrng), cbattr$argvar))

  # Get right function
  fun <- paste0(constr, "Constr")

  #----- Constraint matrix

  if(dim == "var"){

    # The constraint is on var, call directly the method
    Cvar <- do.call(fun, c(list(x = varbasis), pars))

    # For lags, adjust with the rang on other dimension
    odrng <- chkrng(odrng, lagrng, msg = FALSE)
    lagseq <- seq(max(odrng[1], lagrng[1]), min(odrng[2], lagrng[2]), by = 1)
    lagbasis <- do.call(dlnm::onebasis, c(list(x = lagseq), cbattr$arglag))
    Clag <- shapeConstr(lagbasis, shape = "pos", range = odrng)
    Clag$lb <- Clag$ub <- rep(1, NROW(Clag$Cmat))

    # Possibility to compute overall when dim = "var"
    # If so, return directly the result here
    if (overall){
      M <- diag(attr(x, "df")[1]) %x% t(colSums(lagbasis))
      ovCmat <- Cvar$Cmat %*% M
      cmlist <- utils::modifyList(Cvar, list(Cmat = ovCmat))
      return(cmlist)
    }

  # If not var, then do the main constraint on lag dimensions
  } else {

    # Lag part is called directly
    lagbasis <- do.call(dlnm::onebasis, c(list(x = lagrng), cbattr$arglag))
    Clag <- do.call(fun, c(list(x = lagbasis), pars))

    # Var part is adjuste by odrng
    odrng <- chkrng(odrng, varrng, msg = FALSE)
    Cvar <- shapeConstr(varbasis, shape = "pos", range = odrng)
    Cvar$lb <- Cvar$ub <- rep(1, NROW(Cvar$Cmat))
  }

  # Put everything together and return
  Map("%x%", Cvar, Clag)

}
