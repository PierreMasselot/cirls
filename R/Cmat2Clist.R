################################################################################
#
# Create/check list of constraint matrix + bounds from initial list
#
################################################################################

Cmat2Clist <- function(cm, label, nc) {

  # extract objects
  Cmat <- cm$Cmat
  if(!is.null(Cmat)) Cmat <- as.matrix(Cmat)
  lb <- cm$lb
  ub <- cm$ub

  # define nrow
  nr <- ifelse(is.null(Cmat), nc, nrow(Cmat))

  # Cmat: default to simple bound constraints
  if(is.null(Cmat)) Cmat <- diag(nc)
  if(ncol(Cmat) != nc)
    "Constraint matrix for `%s` inconsistent with model frame" |>
    sprintf(label) |> stop()

#########
# NOTE: DEFAULT lb SET TO 0, NOT -Inf
#########

  # bounds
  if(is.null(lb)) {
    # Warning can be a bit annoying. Add a switch
    "No `lb` found for %s. Setting it to 0" |>
    sprintf(label) |> warning()
    lb <- 0
  }
  lb <- rep_len(lb, nr)
  if(is.null(ub)) {
    "No `ub` found for %s. Setting it to Inf" |>
    sprintf(label) |> warning()
    ub <- Inf
  }
  ub <- rep_len(ub, nr)
  if(any(lb > ub))
    "lb greater than ub for `%s`" |> sprintf(label) |> stop()

  # Return the list
  list(Cmat = Cmat, lb = lb, ub = ub)
}
