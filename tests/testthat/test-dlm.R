################################################################################
#
# Tests for constrained DLM
#
################################################################################

library(dlnm)

#-----------------------------
# Simulate data
#-----------------------------

# Some parameters
maxlag <- 21
flag <- function(l) .5 + sin(l / pi)

#----- Simulate

# Use actual data for X
X <- london$tmean
X <- (X - min(X, na.rm = T)) / diff(range(X, na.rm = T))

# Compute cumulated effect
seqlag <- 0:maxlag
Q <- sapply(seqlag, function(l) c(rep(NA, l), X[1:(length(X) - l)]))
cumeff <- Q %*% (flag(0:maxlag) / 5)

# Simulate response
set.seed(24)
suppressWarnings(y <- rpois(length(X), exp((log(20) + cumeff))))

#----- Prepare basis

# A DLM with ns
dlb <- crossbasis(X, lag = maxlag, argvar = list(fun = "lin"),
    arglag = list(fun = "ns", df = 10))

# Test unconstrained
udlm <- glm(y ~ dlb, family = "quasipoisson")
cp <- crosspred(dlb, udlm)

# Plot
plot(cp, ptype = "slices", var = 1)

#-----------------------------
# Test shape constraints
#-----------------------------

test_that("shape-constrained DLM works", {

# Positive
posdlm <- glm(y ~ dlb, family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(dlb, dim = "lag", shape = "pos"))
cp <- crosspred(dlb, posdlm, at = 1)
expect_true(all(cp$matfit > -sqrt(.Machine$double.eps)))

# Decreasing
decdlm <- glm(y ~ dlb, family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(dlb, dim = "lag", shape = "dec"))
cp <- crosspred(dlb, decdlm, at = 1)
expect_true(all(diff(cp$matfit) > -sqrt(.Machine$double.eps)))

# Positive decreasing
pddlm <- glm(y ~ dlb, family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(dlb, dim = "lag", shape = c("pos", "dec")))
cp <- crosspred(dlb, pddlm, at = 1)
expect_true(all(cp$matfit > -sqrt(.Machine$double.eps)))
expect_true(all(diff(cp$matfit) > -sqrt(.Machine$double.eps)))
})

#-----------------------------
# Bound constraints
#-----------------------------

test_that("Bound-constrained DLM works", {

  # Simple zero bound
  bndlm <- glm(y ~ dlb, family = "quasipoisson", method = "cirls.fit",
    constr = ~ bound(dlb, dim = "lag"))
  cp <- crosspred(dlb, bndlm, at = 1)
  expect_equal(cp$matfit[length(cp$matfit)], 0, ignore_attr = TRUE,
    tolerance = 10e-8)

  # Bound and decreasing
  suppressWarnings({
  bndecdlm <- glm(y ~ dlb, family = "quasipoisson", method = "cirls.fit",
    constr = ~ bound(dlb, dim = "lag") + shape(dlb, shape = "dec", dim = "lag"))
  })
  cp <- crosspred(dlb, bndecdlm, at = 1)
  expect_equal(cp$matfit[length(cp$matfit)], 0, ignore_attr = TRUE,
    tolerance = 10e-8)
  expect_true(all(diff(cp$matfit) > -sqrt(.Machine$double.eps)))

  # Bound and positive
  suppressWarnings({
  bnposdlm <- glm(y ~ dlb, family = "quasipoisson", method = "cirls.fit",
    constr = ~ bound(dlb, dim = "lag") + shape(dlb, shape = "pos", dim = "lag"))
  })
  cp <- crosspred(dlb, bnposdlm, at = 1)
  expect_equal(cp$matfit[length(cp$matfit)], 0, ignore_attr = TRUE,
    tolerance = 10e-8)
  expect_true(all(cp$matfit > -sqrt(.Machine$double.eps)))
})

