################################################################################
#
# Create list of constraint matrices from formula
#
################################################################################

constr2Clist <- function(constr, mf = NULL){

  #----- Check terms and data

  # Coerce constr
  constr <- stats::as.formula(constr)

  # Check data
  mf <- mf %||% environment(constr)

  # Terms and dimensions
  term_labs <- names(mf)

  # Add intercept to capture a potential constraint
  term_labs <- c("(Intercept)", term_labs)
  mf[["(Intercept)"]] <- rep(1, NROW(mf))

  # Extract terms for the constr formula
  ct <- attr(stats::terms(constr), "variables")

  # Transform functions if necessary and check if exists
  dropfuns <- rep(FALSE, length(ct))
  for (i in seq(2, length(ct))) {
    ctcstr <- paste0(ct[[i]][[1]], "Constr")
    if (exists(ctcstr)) {
      ct[[i]][[1]] <- str2lang(ctcstr)
    } else {
      dropfuns[i] <- !exists(ct[[i]][[1]])
    }
  }

  # Check variables
  cvars <- lapply(ct[-1], all.vars)
  dropvars <- c(FALSE, sapply(cvars, function(v) any(!v %in% term_labs)))

  # Drop terms
  dropterms <- dropvars | dropfuns
  if (any(dropterms)){
    warning(sprintf("unknown terms in constr have been dropped: %s",
      paste(attr(stats::terms(constr), "variables")[dropterms], collapse = ", ")))
    ct <- ct[!dropterms]
  }
  if (length(ct) == 0) warning("Empty constr list")

  #----- Extract constraint matrices

  # Evaluate the expressions
  Clist <- eval(ct, mf)

  # Add names of the terms
  for (i in seq_along(Clist)) attr(Clist[[i]], "vars") <- cvars[[i]]
  names(Clist) <- attr(stats::terms(constr), "term.labels")[!dropterms[-1]]

  # Return the list
  Clist
}
