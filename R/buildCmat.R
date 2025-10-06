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
buildCmat <- function(mf, assign = NULL, constr = NULL, Cmat = NULL, lb = NULL,
  ub = NULL) {

  # check mf
  if (is.null(attr(mf, "terms"))){
    mf <- stats::model.frame(stats::reformulate(names(mf)), mf)
  }

  # Extract terms in model and check assign
  mt <- attr(mf, "terms")
  termlab <- attr(mt, "term.labels")
  assign <- assign %||% assign(stats::model.matrix(mt, mf))

  # define intercept and first factor (if any), used later
  intercept <- attr(mt, "intercept")
  firstF <- sapply(termlab, function(nm)
    is.factor(mf[[nm]]) || is.logical(mf[[nm]]) || is.character(mf[[nm]]))
  if(sum(firstF) > 1) firstF[which(firstF)[-1]] <- FALSE

  #----- Return empty Cmat/lb/ub if no constraints provided

  # Check if any constraints are passed
  if (all(sapply(list(Cmat, lb, ub, constr), is.null))) {
    Cmat <- matrix(nrow = 0, ncol = length(assign))
    lb <- ub <- numeric(0)
    cmempty <- list(Cmat = Cmat, lb = lb, ub = ub)
    return(cmempty)
  }

  #----- Check and return Cmat/lb/ub if provided in full

  # if Cmat/lb/ub are numeric, they represent the full constraints
  if (any(sapply(list(Cmat, lb, ub), is.numeric))) {
    cmfull <- Cmat2Clist(list(Cmat, lb, ub), label="Cmat", nc=length(assign))
    return(cmfull)
  }

  #----- Prepare constraints related to Cmat/lb/ub

  # Identify the terms and check if in model formula
  cmterms <- unique(c(names(Cmat), names(lb), names(ub)))
  if(any(ind <- !cmterms %in% termlab))
    stop(sprintf("term(s) in Cmat/lb/ub not in model formula: %s",
      paste(cmterms[ind], collapse = ", ")))

  # create the Cmat/lb/ub list for each term
  cmlist <- lapply(cmterms, function(nm)
    list(Cmat=Cmat[[nm]], lb=lb[[nm]], ub=ub[[nm]]))
  names(cmlist) <- cmterms

  # extract ncols for terms
  cmnc <- sapply(match(cmterms, termlab), function(x) sum(assign==x))

  # Check and complete the list
  cmlist <- mapply(Cmat2Clist, cm = cmlist, label = cmterms, nc = cmnc,
    SIMPLIFY = FALSE)

  #----- Prepare constraints related to constr

  # Coerce constr, extract expressions
  constr <- if(!is.null(constr)) stats::as.formula(constr)
  csvars <- if(!is.null(constr))
    attr(stats::terms(constr), "variables") else list()

  # Identify the terms (optionally multiple) and check if in model formula
  csterms <- all.vars(csvars, unique=FALSE)
  if(any(ind <- !csterms %in% termlab))
    stop(sprintf("term(s) in constr not in model formula: %s",
      paste(csterms[ind], collapse = ", ")))

  # Intercept indicator
  intind <- if(!intercept && any(firstF)) csterms %in% termlab[firstF] else
    rep_len(F, length(csterms))

  # create the Cmat/lb/ub list for each term
  cslist <- mapply(constr2Clist, var = csvars[-1], label = csterms,
    int = intind, MoreArgs = list(mf = mf), SIMPLIFY = FALSE)
  names(cslist) <- csterms

  #----- Create the full Cmat/lb/ub

  # Initialise objects
  alllist <- c(cmlist, cslist)
  allterms <- c(cmterms, csterms)
  nr <- sapply(alllist, function(x) nrow(x$Cmat))
  Cmat <- matrix(0, sum(nr), length(assign))
  lb <- ub <- rep_len(0, sum(nr))

  # Fill
  cnr <- c(0, cumsum(nr))
  for(i in seq(alllist)) {
    indc <- assign == (which(termlab==allterms[i]))
    indr <- seq(cnr[i]+1, cnr[i+1])
    Cmat[indr, indc] <- alllist[[i]]$Cmat
    lb[indr] <- alllist[[i]]$lb
    ub[indr] <- alllist[[i]]$ub
  }

  #----- Final checks and return

  # Check Cmat is irreducible
  if (nrow(Cmat) > 1){
    chkc <- checkCmat(Cmat)
    if (sum(chkc$redundant) > 0){
      Cmat <- Cmat[!chkc$redundant,,drop = F]
      lb <- lb[!chkc$redundant]
      ub <- ub[!chkc$redundant]
      warning(paste0("Redundant constraints removed from Cmat: ",
        paste(which(chkc$redundant), collapse = ", ")))
    }
    if (sum(chkc$equality) > 0){
      warning(paste0("Underlying equality constraints: ",
        paste(which(chkc$equality), collapse = ", "), ". ",
        "Consider using lb and ub to set equality constraints instead."))
    }
  }

  # Return
  list(Cmat = Cmat, lb = lb, ub = ub)
}
