################################################################################
#
# Main function to build Cmat
#
################################################################################

#' Build a constraint matrix
#'
#' @description
#' Internal function building a full constraint matrix from a list of constraint matrices and/or a formula providing specific constraints. In-depth technical details are provided below.
#'
#' @param mf A [model.frame][stats::model.frame()] or a list of variables. Defines the model for which the constraint matrix is built.
#' @param assign A vector of indices mapping columns of the design matrix to model terms.
#' @param constr A formula specifying constraints.
#' @param Cmat Either a matrix or a named list of constraint matrices. In the latter case, names should correspond to terms in `mf` (dee setails).
#' @param lb,ub Vector or named list of vectors containing constraint bounds. In the latter case, names should correspond to terms in `mf` (dee setails).
#'
#' @details
#' This function is called internally by [cirls.fit][cirls.fit()] and is not meant to be used directly by the user. It prepares the full `Cmat`, `lb` and `ub` for the model, providing a way to specify constraints without having to build a full constraint matrix beforehand. It uses the model frame in `mf` to match specific constraints to the right columns in the design matrix.
#'
#' This function also checks that any constraint matrix provided to [cirls.fit][cirls.fit()] is irreducible. See [checkCmat][checkCmat()] for details.
#'
#' There are three ways to specify constraints in `cirls`. Each one is detailed in a subsetction below. Note that options 1 and 2 can be used simultaneously.
#' 1. The `constr` formula provides a simple interface for many commonly encountered constraints. This is the recommended way for new users.
#' 2. A named list of constraint matrices and bounds can be provided to `Cmat`/`lb`/`ub` arguments. This is useful for constraints that are not (yet) implemented for the `constr` formula.
#' 3. A fully specified `Cmat`/`lb`/`ub` can be provided directly.
#'
#' ## The `constr` formula
#'
#' The easiest way to specify constraints is as a formula of the form `constr = ~ cons(x, ...)` where `cons` represents a constraint to be applied to term `x` of the `cirls` model and `...` represents additional arguments that depend on the `cons` function. Several constraints can be applied to the same or to different terms, for example `constr = ~ cons1(x) + cons1(y) + cons2(x)`. All terms appearing in `constr` can either be vectors or matrices and should be found in the model.frame `mf`, otherwise the constraint is dropped with a warning.
#'
#' Internally, `buildCmat` will look for a function called `consConstr` that returns a list of `Cmat`, `lb` and `ub` built specifically for the term provided to `cons()` in the formula. This provides a simple way to add new constraints to the interface, by creating a function with the `Constr` suffix, taking a term as an input and outputting a list with `Cmat`, `lb` and `ub`.
#'
#' The list of available `cons` functions can be found on the main help page of the [cirls][cirls-package] package. Each function has its own help page detailing the implemented constraint with available parameter.
#'
#' ## Term-specific `Cmat`/`lb`/`ub`
#'
#' Constraints for specific terms of the model can be provided as named lists to one or several of arguments `Cmat`,`lb` and `ub` which represent the constraint matrix, lower bound and upper bound respectively. This will take the form `Cmat = list(x = cm1, y = cm2)` where `x` and `y` are terms in the model.frame `mf`, while `cm1` and `cm2` are constraint matrices (or bound vectors for `lb` and `ub`) *that should be consistent with the dimensions of their respective terms*.
#'
#' When terms are found in `lb` and/or `ub`, but not in `Cmat`, it is assumed that *simple bound constraints are to be applied to the terms*. In that case, `Cmat` will be internally created as a simple identity matrix matching the dimensions of the terms in question.
#'
#' Names in `Cmat`/`lb`/`ub` can include several terms, which should be separated by a semicolon `;`, for instance `Cmat = list("x;y" = cm)`. This allows specifying constraints that span several terms in the model.
#'
#' ## Fully specified `Cmat`/`lb`/`ub`
#'
#' If one of `Cmat`/`lb`/`ub` is neither `NULL` nor a list, it is assumed it is provided fully specified, i.e. should be returned as is. In that case, `constr` and any list provided to other arguments are ignored. When not all of `Cmat`/`lb`/`ub` are fully specified, the other ones will be filled with default values that match the dimension of the model matrix, i.e. an identity matrix for `Cmat`, a vector of zeros for `lb` and a vector of `Inf` for `ub`.
#'
#' ## Unconstrained model
#'
#' When all of `constr`/`Cmat`/`lb`/`ub` are `NULL`, a list of empty `Cmat`/`lb`/`ub` is returned and an unconstrained model is fit.
#'
#' @returns A list with elements `Cmat`, `lb` and `ub` containing the fully specified constraint matrix, lower and upper bounds for the model specified in argument `mf`. `Cmat` additionally include an attribute called `terms` which maps constraints represented in the matrix to individual terms in the model.
#'
#' @seealso The main [help page][cirls-package] for the list of `cons` functions implemented and examples. [checkCmat][checkCmat()] for details on irreducibility.
#'
#' @example inst/examples/ex_london_buildCmat.R
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
  assign <- assign %||% attr(stats::model.matrix(mt, mf), "assign")

  # define intercept and first factor (if any), used later
  intercept <- attr(mt, "intercept")
  firstF <- sapply(termlab, function(nm)
    is.factor(mf[[nm]]) || is.logical(mf[[nm]]) || is.character(mf[[nm]]))
  if(sum(firstF) > 1) firstF[which(firstF)[-1]] <- FALSE

  #----- Return empty Cmat/lb/ub if no constraints provided

  # Check if any constraints are passed
  if (all(sapply(list(Cmat, lb, ub, constr), is.null))) {
    warning("No constraint provided")
    Cmat <- matrix(nrow = 0, ncol = length(assign))
    lb <- ub <- numeric(0)
    cmempty <- list(Cmat = Cmat, lb = lb, ub = ub)
    return(cmempty)
  }

  #----- Check and return Cmat/lb/ub if provided in full

  # if Cmat/lb/ub are numeric, they represent the full constraints
  if (any(sapply(list(Cmat, lb, ub), is.numeric))) {
    cmfull <- Cmat2Clist(list(Cmat = Cmat, lb = lb, ub = ub), label = "Cmat",
      nc = length(assign))
    return(cmfull)
  }

  #----- Prepare constraints related to Cmat/lb/ub

  # Identify the terms and check if in model formula
  cmlabs <- unique(c(names(Cmat), names(lb), names(ub)))
  cmterms <- strsplit(as.character(cmlabs), ";")
  ind <- !sapply(cmterms, function(x) all(x %in% termlab))
  if(any(ind)){
    warning(sprintf(
      "dropped term(s) in Cmat/lb/ub that contain variables not in model formula: %s",
      paste(cmlabs[ind], collapse = ", ")))
    cmlabs <- cmlabs[!ind]
    cmterms <- cmterms[!ind]
  }

  # create the Cmat/lb/ub list for each term
  cmlist <- lapply(cmlabs, function(nm)
    list(Cmat = Cmat[[nm]], lb = lb[[nm]], ub = ub[[nm]]))
  names(cmlist) <- cmlabs

  # extract ncols for terms
  cmnc <- sapply(cmterms, function(x) sum(sapply(match(x, termlab),
    function(y) sum(assign==y))))

  # Check and complete the list
  cmlist <- mapply(Cmat2Clist, cm = cmlist, label = cmlabs, nc = cmnc,
    SIMPLIFY = FALSE)

  #----- Prepare constraints related to constr

  # Coerce constr, extract expressions
  constr <- if(!is.null(constr)) stats::as.formula(constr)
  csvars <- if(!is.null(constr))
    attr(stats::terms(constr), "variables") else list()
  cslabs <- sapply(csvars, deparse)[-1]

  # Identify the terms (optionally multiple) and check if in model formula
  csterms <- lapply(csvars[-1], all.vars)
  ind <- !sapply(csterms, function(x) all(x %in% termlab))
  if(any(ind)){
    warning(sprintf(
      "dropped term(s) in constr that contain variables not in model formula: %s",
      paste(cslabs[ind], collapse = ", ")))
    csvars <- csvars[-(which(ind) + 1)]
    csterms <- csterms[!ind]
    cslabs <- cslabs[!ind]
  }

  # Intercept indicator
  intind <- if(!intercept && any(firstF))
    sapply(csterms, "%in%", termlab[firstF]) else
    rep_len(F, length(csterms))

  # create the Cmat/lb/ub list for each term
  cslist <- mapply(constr2Clist, var = csvars[-1], label = cslabs,
    int = intind, MoreArgs = list(mf = mf), SIMPLIFY = FALSE)
  names(cslist) <- cslabs

  #----- Create the full Cmat/lb/ub

  # Initialise objects
  alllist <- c(cmlist, cslist)
  allterms <- c(cmterms, csterms)
  if (length(alllist) == 0) warning("No valid constraint provided")
  nr <- sapply(alllist, function(x) nrow(x$Cmat)) |> unlist()
  Cmat <- matrix(0, sum(nr), length(assign))
  lb <- ub <- rep_len(0, sum(nr))
  rownames(Cmat) <- names(lb) <- names(ub) <- rep(names(alllist), nr)

  # Fill
  cnr <- c(0, cumsum(nr))
  for(i in seq(alllist)) {
    indc <- lapply(match(allterms[[i]], termlab),
      function(x) which(assign == x)) |> unlist()
    indr <- seq(cnr[i] + 1, cnr[i + 1])
    Cmat[indr, indc] <- alllist[[i]]$Cmat
    lb[indr] <- alllist[[i]]$lb
    ub[indr] <- alllist[[i]]$ub
  }

  #----- Final checks and return

  # Reduce Cmat
  redCm <- checkCmat(Cmat = Cmat, lb = lb, ub = ub, reduce = TRUE, warn = TRUE)

  # Return
  redCm[c("Cmat", "lb", "ub")]
}
