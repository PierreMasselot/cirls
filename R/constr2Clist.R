################################################################################
#
# Create list of constraint matrices from formula
#
################################################################################

constr2Clist <- function(constr, mf = NULL){

  #----- Extract data and terms

  # Coerce constr
  constr <- stats::as.formula(constr)

  # Check data
  mf <- mf %||% environment(constr)

  # Terms and dimensions
  term_labs <- names(mf)

  # Add intercept to capture a potential constraint
  intercept <- attr(attr(mf, "terms"), "intercept") %||% FALSE
  if (intercept){
    term_labs <- c("(Intercept)", term_labs)
    mf[["(Intercept)"]] <- rep(1, NROW(mf))
  }

  # Extract terms for the constr formula
  ct <- attr(stats::terms(constr), "variables")

  # Check if there are factors
  isF <- sapply(mf,
    function(x) is.factor(x) || is.logical(x) || is.character(x))
  firstF <- term_labs[which(isF)[1]]

  #----- Check data and functions can be matched to existing objects

  # Initialise a vector to drop terms that cannot be interpreted
  # And another one to get the variables involved in each term
  dropterms <- rep(FALSE, length(ct))
  cvars <- rep("", length(ct) - 1)

  # Loop across the terms
  for (i in seq(2, length(ct))) {

    # Check if the constraint function exists
    ctcstr <- paste0(ct[[i]][[1]], "Constr")
    if (exists(ctcstr)) {
      ct[[i]][[1]] <- str2lang(ctcstr)
    } else {
      dropterms[i] <- !exists(ct[[i]][[1]])
    }

    # Check if variables can be found in the model frame
    cvars <- all.vars(ct[[i]])
    # NB: logicals 'T' and 'F' are treated as variables in R
    if (!all(cvars %in% c(term_labs, "T", "F"))){
      dropterms[i] <- TRUE
    } else {
      # If the terms include the first factor, add 'intercept' argument
      # Because factor will be coded as dummy regardless of contrasts
      # If such an argument doesn't exist, it will be ignored
      if (any(cvars %in% firstF) && !intercept){
        ct[[i]]$intercept <- TRUE
      }
    }
  }

  # Drop terms that have been flagged with a warning
  if (any(dropterms)){
    warning(paste0("unknown terms in 'constr' have been dropped: ",
      paste(attr(stats::terms(constr), "variables")[dropterms],
        collapse = ", "),
      "\nIt is recommended to check the resulting 'Cmat'."))
    ct <- ct[!dropterms]
  }
  if (length(ct) == 0) warning("Empty 'constr' list")

  #----- Extract constraint matrices

  # Evaluate the expressions
  Clist <- eval(ct, mf)

  # Add names of the terms
  for (i in seq_along(Clist)) attr(Clist[[i]], "vars") <- all.vars(ct[[i + 1]])
  names(Clist) <- attr(stats::terms(constr), "term.labels")[!dropterms[-1]]

  # Return the list
  Clist
}
