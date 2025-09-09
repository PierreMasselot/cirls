################################################################################
#
# Main function to build Cmat
#
################################################################################

#### Perhaps keep a help page here for technical details. Attach the usage to a clearer help page

#' Build a constraint matrix
#'
#' @description
#' Function building a full constraint matrix from a list of constraint matrices and/or a formula providing specific constraints. Mainly used internally by [cirls.fit][cirls.fit()].
#'
#' @param mf A [model.frame][stats::model.frame()] or a list of variables.
#' @param constr A formula specifying constraints.
#' @param Cmat A named list of constraint matrices where names should be found among the terms in `mf`.
#' @param lb,ub Vector or list of vectors containing constraint bounds. If a vector, is used as default bounds for terms with no specified bounds. If a named list, is matched to `Cmat` to provide corresponding bounds.
#'
#' @details
#' This function is called internally by `cirls.fit` whenever `Cmat` is not a matrix, and provides a way to specify constraints without having to build a full constraint matrix beforehand. It uses the model frame in `mf` to match specific constraints to the right columns in the design matrix.
#'
#' The argument `constr` provides a simple way to specify potentially complex constraints. It is a formula of the form `~ shape(x, ...)` where `shape` specifies the constraint and `x` the term in `mf` to which it applies. Internally, the formula will look for a function named `shapeConstr` to be called on variable `x` (which allows for several columns). The `...` represent potential additional arguments for the `shapeConstr` function. For the list of available constraints and how to create new ones, see **upcoming**.
#'
#' The argument `Cmat` is used to provide a named list of constraint matrices, where names should correspond to terms in `mf`. This allows providing custom constraint matrices to specific terms that wouldn't be available through `constr`. Names in `Cmat` can include several terms, which should be separated by a `;`, for instance `x1;x2`. Although not mandatory, elements in `Cmat` can have attributes `lb`, `ub` and `vars` to provide lower and upper bounds, and term names, respectively.
#'
#' `lb` and `ub` are meant to be used in conjunction with `Cmat`. If a simple value or vector, they will be used as default values for elements in `Cmat` for which no bounds is specified in its attributes. If lists, they provide bounds for constraint matrices in `Cmat`. In this case, all the names in `Cmat` should be found in `lb` and `ub`.
#'
#' Note that both `constr` and `Cmat` can be used at the same time, and neither is mandatory. If both are `NULL`, an empty constraint matrix will be returned.
#'
#' @returns A list with containing elements `Cmat`, `lb` and `ub` containing the full constraint matrix, lower and upper bounds for the model specified in argument `mf`. `Cmat` additionally include an attribute called `terms` which maps constraints represented in the matrix to individual terms in the model.
#'
#' @examples
#' ####### Upcoming
#'
#' @export
buildCmat <- function(mf, constr = NULL, Cmat = NULL, lb = 0, ub = Inf){

  # If not a model frame, try to coerce to one (assuming a simple model)
  if (is.null(attr(mf, "terms"))){
    mf <- stats::model.frame(stats::reformulate(names(mf)), mf)
  }

  # Extract terms in model
  mt <- attr(mf, "terms")
  term_labs <- attr(mt, "term.labels")
  tcols <- stats::model.matrix(mt, mf) |> attr("assign")

  # Add info about intercept
  if (attr(mt, "intercept")){
    term_labs <- c("(Intercept)", term_labs)
  }

  #----- Prepare constraints

  # Go through Cmat to attributing them to terms
  dropterms <- rep(FALSE, length(Cmat))
  for (i in seq_along(Cmat)){
    nm <- names(Cmat)[i]
    attr(Cmat[[i]], "vars") <- unlist(strsplit(attr(Cmat[[i]], "vars") %||% nm,
      ";"))
    dropterms[i] <- any(!attr(Cmat[[i]], "vars") %in% term_labs)
  }

  # Drop terms that weren't found (with warning)
  if (any(dropterms)){
    warning(sprintf("Unknown terms in Cmat have been dropped: %s",
      paste(names(Cmat)[dropterms], collapse = ", ")))
    Cmat <- Cmat[!dropterms]
  }

  # Add constraints in constr formula
  if (!is.null(constr)){
    Cmat <- c(Cmat, constr2Clist(constr, mf))
  }

  # Check there are constraints left
  if (length(Cmat) == 0) {
    warning("No valid constraint found. Fitting unconstrained model")
    return(list(Cmat = matrix(nrow = 0, ncol = length(tcols)),
      lb = numeric(0), ub = numeric(0)))
  }

  #----- Build full constraint matrix

  # Prepare lbs and ubs
  if(!is.list(lb)){
    lb <- rep(list(lb), length(Cmat))
    names(lb) <- names(Cmat)
  }
  if(!is.list(ub)){
    ub <- rep(list(ub), length(Cmat))
    names(ub) <- names(Cmat)
  }

  # Go through each term in Clist to expand constraint matrices
  cml <- lbl <- ubl <- vector("list", length(Cmat))
  names(cml) <- names(lbl) <- names(ubl) <- names(Cmat)
  for (i in seq_along(Cmat)) {

    # Initialise a matrix with 0 values
    nc <- NROW(Cmat[[i]])
    cml[[i]] <- matrix(0, nc, length(tcols))

    # Find the columns to replace
    cmterms <- match(attr(Cmat[[i]], "vars"), term_labs)
    cols <- tcols %in% (cmterms - attr(mt, "intercept"))

    # Check if the dimensions match and replace
    if (sum(cols) == NCOL(Cmat[[i]])) cml[[i]][, cols] <- Cmat[[i]] else stop(
      sprintf("Constraint matrix for term `%s` inconsistent with model frame",
        names(Cmat)[i]))

    # Match bounds
    lbl[[i]] <- attr(Cmat[[i]], "lb") %||% lb[[names(Cmat)[i]]] %||% {
      warning(sprintf("No `lb` found for constraint %s. Setting it to -Inf",
        names(Cmat)[i]))
      -Inf
    } |> rep_len(nc)
    ubl[[i]] <- attr(Cmat[[i]], "ub") %||% ub[[names(Cmat)[i]]] %||% {
      warning(sprintf("No `ub` found for constraint %s. Setting it to Inf",
        names(Cmat)[i]))
      Inf
    } |> rep_len(nc)
  }

  # Append terms together
  Cmat <- do.call("rbind", cml)
  lb <- unlist(lbl)
  ub <- unlist(ubl)

  #----- Tidy and return

  # Check Cmat is irreducible
  if (nrow(Cmat) > 1){
    chkc <- checkCmat(Cmat)
    if (length(chkc$redundant) > 0){
      Cmat <- Cmat[-chkc$redundant,,drop = F]
      lb <- lb[-chkc$redundant]
      ub <- ub[-chkc$redundant]
      warning(paste0("Redundant constraints removed from Cmat: ",
        paste(chkc$redundant, collapse = ", ")))
    }
    if (length(chkc$equality) > 0){
      warning(paste0("Underlying equality constraints: ",
        paste(chkc$equality, collapse = ", "), ". ",
        "Consider using lb and ub to set equality constraints instead."))
    }
  }

  # Add colnames and add terms as attribute
  colnames(Cmat) <- colnames(stats::model.matrix(mt, mf))
  rownames(Cmat) <- names(lb)

  # Attributes constraints to lines
  nconsterm <- sapply(cml, NROW)
  end <- cumsum(nconsterm)
  start <- end - nconsterm + 1
  conslist <- Map(seq, start, end)
  attributes(Cmat)$terms <- conslist

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
