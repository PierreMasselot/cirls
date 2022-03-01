# Transform list of constraint matrices as one big constraint matrix
clist2cmat <- function(Clist, mt, x)
{
  Clist <- lapply(Clist, as.matrix)
  #----- Extract number of variables in each term
  # Extract variable names
  nterms <- length(attr(mt, "term.labels")) + attr(mt, "intercept")
  term_labs <- attr(mt, "term.labels")
  # Add intercept if present
  if (attr(mt, "intercept")) term_labs <- c("(Intercept)", term_labs)
  ncolterm <- table(attr(x, "assign"))
  #----- Number of constraint for each term
  # Check the names of Clist
  if (any(!names(Clist) %in% term_labs))
    stop("unknown terms in Cmat")
  # 0 by default
  nconsterm <- rep(0, nterms)
  # Fill for terms with provided constrain matrix
  whichcons <- match(names(Clist), term_labs)
  nconsterm[whichcons] <- sapply(Clist, NROW)
  # Check that number of columns match
  diffncol <- sapply(Clist, ncol) != ncolterm[whichcons]
  if(any(diffncol)){
    stop(
      sprintf("constraint matrix of term(s) %s has the wrong number of columns",
        paste(term_labs[diffncol], sep = ", "))
    )
  }
  #----- Fill the big constraint matrix
  # Get where to replace matrix
  cdim <- cbind(nconsterm, ncolterm)
  end <- apply(cdim,2,cumsum)
  start <- end - cdim + 1
  matind <- array(seq(prod(colSums(cdim))),colSums(cdim))
  ind <- unlist(lapply(whichcons,function(i)
    matind[start[i,1]:end[i,1],start[i,2]:end[i,2]]))
  # Initialize with zeros and replace
  cmat <- matrix(0, sum(nconsterm), sum(ncolterm))
  cmat[ind] <- unlist(Clist)
  # Return
  cmat
}
