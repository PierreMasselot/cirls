################################################################################
#
# Check if constraint matrix is reducible
# Returns constraints that can be removed safely
#
################################################################################

check_cmat <- function(Cmat){
  tCmat <- t(Cmat)
  ncons <- ncol(tCmat)
  redundant <- equality <- rep(FALSE, ncons)
  for (i in seq_len(ncons)){
    y <- tCmat[,i]
    x <- tCmat[,-c(i, which(redundant)),drop = F]
    # Check redundancy
    res <- coneproj::coneB(y, x)
    redundant[i] <- isTRUE(all.equal(y, drop(res$yhat)))
    if (!redundant[i]){
      # Check underlying equality constraint
      reseq <- limSolve::nnls(x, -y)
      equality[i] <- isTRUE(all.equal(-y, drop(x %*% reseq$X)))
    }
  }
  # Return indices
  list(redundant = which(redundant), equality = which(equality))
}
