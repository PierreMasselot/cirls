################################################################################
#
#  Test main routine
#
################################################################################

#------------------------------
# Simple test
#------------------------------

#----- Parameters

# Number of obs
n <- 1000

# Coefficients
betas <- c(0, 1, 2, -1, 1)
p <- length(betas)

#----- Generate data

# Uniform values between 0 and 1
x <- matrix(rnorm(n * p), n, p)

# Linear predictor
eta <- 5 + x %*% betas

# Simulate responses
ynorm <- eta + rnorm(n, 0, .2)
ypois <- rpois(n, exp(eta))

#----- Different constraint matrices

# Everything positive
cpos <- diag(p)

# Increasing
cinc <- diff(diag(p))

#----- Apply models
normpos <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos))
norminc <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cinc))
poispos <- glm(ypois ~ x, family = "poisson",
  method = cirls.fit, Cmat = list(x = cpos))
poisinc <- glm(ypois ~ x, family = "poisson",
  method = cirls.fit, Cmat = list(x = cinc))

#----- Tests ------

test_that("output contains all new elements", {

  # Contains all new outputs
  expect_true(all(
    c("active.cons", "Cmat", "lb", "ub") %in%
    names(normpos)))
  expect_true(all(
    c("active.cons", "Cmat", "lb", "ub") %in%
    names(norminc)))
  expect_true(all(
    c("active.cons", "Cmat", "lb", "ub") %in%
    names(poispos)))
  expect_true(all(
    c("active.cons", "Cmat", "lb", "ub") %in%
    names(poisinc)))

  # Cmat is correct
  expect_equal(normpos$Cmat, cbind(0, cpos))
  expect_equal(norminc$Cmat, cbind(0, cinc))
  expect_equal(poispos$Cmat, cbind(0, cpos))
  expect_equal(poisinc$Cmat, cbind(0, cinc))

  # bounds are well computed
  expect_equal(normpos$lb, rep(0, p))
  expect_equal(norminc$lb, rep(0, p - 1))
  expect_equal(poispos$lb, rep(0, p))
  expect_equal(poisinc$lb, rep(0, p - 1))
  expect_equal(normpos$ub, rep(Inf, p))
  expect_equal(norminc$ub, rep(Inf, p - 1))
  expect_equal(poispos$ub, rep(Inf, p))
  expect_equal(poisinc$ub, rep(Inf, p - 1))

})

test_that("results respect constraints", {
  # positive
  expect_true(all(coef(normpos)[-1] >= (0 - 1e-6)))
  expect_true(all(coef(poispos)[-1] >= (0 - 1e-6)))

  # increases
  expect_true(all(diff(coef(norminc)[-1]) >= (0 - 1e-6)))
  expect_true(all(diff(coef(poisinc)[-1]) >= (0 - 1e-6)))
})


#----- Test lb and ub ------

# Revert constraints
normpos_rev <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = -cpos),
  lb = -Inf, ub = 0)
norminc_rev <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = -cinc),
  lb = -Inf, ub = 0)

# Check this is equaivalent to previous constraint
test_that("reverting ub and lb is equal", {
  expect_mapequal(coef(normpos), coef(normpos_rev))
  expect_mapequal(coef(norminc), coef(norminc_rev))
})

# Equality constraint
betasum <- sum(betas)
normeq <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = betasum, ub = betasum)
test_that("equality constraint can be passed", {
  expect_equal(sum(coef(normeq)[-1]), betasum)
})

#----- Alternative solvers ------

# quadprog
quadprog_pos <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  qp_solver = "quadprog")
quadprog_inc <- glm(ypois ~ x, method = cirls.fit, Cmat = list(x = cinc),
  qp_solver = "quadprog")
betasum <- sum(betas)
quadprog_eq <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = betasum, ub = betasum, qp_solver = "quadprog")

test_that("quadprog solver is integrated", {
  expect_true(all(coef(quadprog_pos)[-1] >= (0 - 1e-6)))
  expect_true(all(diff(coef(quadprog_inc)[-1]) >= (0 - 1e-6)))
  expect_equal(sum(coef(quadprog_eq)[-1]), betasum)
})

# coneproj
cone_pos <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  qp_solver = "coneproj")
cone_inc <- glm(ypois ~ x, method = cirls.fit, Cmat = list(x = cinc),
  qp_solver = "coneproj")
betasum <- sum(betas)
cone_eq <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = betasum, ub = betasum, qp_solver = "coneproj")

test_that("coneproj solver is integrated", {
  expect_true(all(coef(cone_pos)[-1] >= (0 - 1e-6)))
  expect_true(all(diff(coef(cone_inc)[-1]) >= (0 - 1e-6)))
  expect_equal(sum(coef(cone_eq)[-1]), betasum)
})
