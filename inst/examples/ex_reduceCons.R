###################################################
# Simple example

# Two constraints: one >=0 and one >= 1
Cmat <- as.matrix(rep(1, 2))
lb <- 0:1

# The first one is removed because redundant
reduceCons(Cmat, lb)

###################################################
# Example of underlying equality constraint

# Contraint: Parameters sum is >= 0 and sum is <= 0
cmateq <- rbind(rep(1, 3), rep(-1, 3))

# Transformed into a single equality constraint
reduceCons(cmateq)

###################################################
# An example with shape constraints

library(dlnm)

# Constraints: successive coefficients should increase and be convex
# Intuitively, if the first two coefficients increase,
# then convexity forces the rest to increase which means there is redundancy
basis <- ps(1:10, df = 5)
cmatic <- rbind(
  shapeConstr(basis, shape = "inc")$Cmat, # Increasing
  shapeConstr(basis, shape = "cvx")$Cmat # Convex
)

# Checking indicates that some constraints are redundant
# Returns reduced matrix and a warning
reduceCons(cmatic)

# Compare without removing the redundant constraints
reduceCons(cmatic, redundant = FALSE)

# Note that this is silently done when both "inc" and "cvx" are provided
shapeConstr(basis, shape = c("inc", "cvx"))
