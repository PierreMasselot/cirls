# Transform list of constraint matrices as one big constraint matrix
clist2cmat <- function(Clist){
  # Get the objects from parent environments
  # that cannot be passed through parameters
  mt <- get("mt", envir = parent.frame(3)) 
  nterms <- length(attr(mt, "term.labels")) + attr(mt, "intercept")
  term_labs <- attr(mt, "term.labels")
  if (attr(mt, "intercept")) term_labs <- c("(Intercept)", term_labs)
  ncolterm <- table(attr(get("x", envir = parent.frame(2)), "assign"))
  # Check the names of Clist
  if (any(!names(Clist) %in% term_labs))
    stop("unknown terms in Cmat")
  # Create a list of matrices with each term in the right order
  matlist <- rep(list(NULL), nterms)
  # Add the provided matrices for each term
  matlist[match(names(Clist), term_labs)] <- Clist
  # Create empty matrices for missing terms
  empty <- lengths(matlist) == 0
  matlist[empty] <- lapply(ncolterm[empty], 
    function(n) matrix(nrow = 0, ncol = n))
  # Check that the number of columns matches the terms
  diffncol <- sapply(matlist, ncol) != ncolterm
  if(any(diffncol)){
    stop(
      sprintf("constraint matrix of term(s) %s has the wrong number of columns",
        paste(term_labs[diffncol], sep = ", "))
    )
  }
  # Put all matrices in a big constraint one
  as.matrix(Matrix::bdiag(matlist))
}