################################################################################
#
# Transform list of constraint matrices as one big constraint matrix
#
################################################################################

clist2cmat <- function(Clist, mt, x)
{
  # prepare matrices
  Clist <- lapply(Clist, as.matrix)

  #----- Extract number of variables in each term

  # Extract variable names
  nterms <- length(attr(mt, "term.labels")) + attr(mt, "intercept")
  term_labs <- attr(mt, "term.labels")

  # Add intercept if present
  if (attr(mt, "intercept")) term_labs <- c("(Intercept)", term_labs)
  ncolterm <- table(attr(x, "assign"))

  #----- Number of constraint for each term

  # Check the names of Clist and drop unknown ones
  if (any(!names(Clist) %in% term_labs)){
    warning("unknown terms in Cmat have been dropped")
    Clist <- Clist[names(Clist) %in% term_labs]
  }
  if (length(Clist) == 0) stop("Empty Cmat list")

  # 0 by default
  nconsterm <- rep(0, nterms)

  # Fill for terms with provided constraint matrix
  whichcons <- match(names(Clist), term_labs)
  nconsterm[whichcons] <- sapply(Clist, NROW)

  # Check that number of columns match
  diffncol <- sapply(Clist, ncol) != ncolterm[whichcons]
  if(any(diffncol)){
    stop(
      sprintf("constraint matrix of term(s) %s has the wrong number of columns",
        paste(term_labs[diffncol], collapse = ", "))
    )
  }

  #----- Fill the big constraint matrix

  # Get where to replace matrix
  cdim <- cbind(nconsterm, ncolterm)
  end <- cbind(cumsum(nconsterm), cumsum(ncolterm))
  start <- end - cdim + 1
  matind <- array(seq(prod(colSums(cdim))), colSums(cdim))
  ind <- unlist(lapply(whichcons, function(i)
    matind[start[i,1]:end[i,1], start[i,2]:end[i,2]]))

  # Initialize with zeros and replace
  Cmat <- matrix(0, sum(nconsterm), sum(ncolterm))
  Cmat[ind] <- unlist(Clist)

  # Add colnames and add terms as attribute
  colnames(Cmat) <- dimnames(x)[[2L]]
  conslist <- lapply(whichcons, function(i) start[i,1]:end[i,1])
  names(conslist) <- names(Clist)
  attributes(Cmat)$terms <- conslist

  # Return
  Cmat
}
