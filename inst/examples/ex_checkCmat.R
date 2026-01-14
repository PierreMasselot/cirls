###################################################
# Example of reducible matrix

# Constraints: successive coefficients should increase and be convex
# Intuitively, if the first two coefficients increase,
# then convexity forces the rest to increase which means there is redundancy
p <- 5
cmatic <- rbind(
  shapeConstr(matrix(NA, 0, p), shape = "inc")$Cmat, # Increasing
  shapeConstr(matrix(NA, 0, p), shape = "cvx")$Cmat # Convex
)

# Checking indicates that some constraints are redundant.
# Returns reduced matrix and a warning
checkCmat(cmatic)

# Note that this is silently done when both "inc" and "cvx" are provided
shapeConstr(matrix(NA, 0, p), shape = c("inc", "cvx"))$Cmat

###################################################
# Example of irreducible matrix

# Constraints: coefficients form an S-shape
p <- 4
cmats <- rbind(
  diag(p)[1,], # positive
  diff(diag(p))[c(1, p - 1),], # Increasing at both end
  diff(diag(p), diff = 2)[1:(p/2 - 1),], # First half convex
  -diff(diag(p), diff = 2)[(p/2):(p-2),] # second half concave
)

# Note, this matrix is not of full row rank
qr(t(cmats))$rank
all.equal(cmats[2,] + cmats[4,] - cmats[5,], cmats[3,])

# However, it is irreducible: all constraints are necessary
checkCmat(cmats)

###################################################
# Example of underlying equality constraint

# Contraint: Parameters sum is >= 0 and sum is <= 0
cmateq <- rbind(rep(1, 3), rep(-1, 3))

# Checking indicates that both constraints imply equality constraint
checkCmat(cmateq)
