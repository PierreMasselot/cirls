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
set.seed(1234)
x <- matrix(rnorm(n * p), n, p)

# Linear predictor
eta <- 1 + x %*% betas

# Simulate responses
set.seed(5678)
ynorm <- eta + rnorm(n, 0, .2)
ypois <- rpois(n, exp(eta))

#----- Different constraint matrices

# Everything positive
cpos <- diag(p)

# Increasing
cinc <- diff(diag(p))

#----- Apply models
normpos <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = 0), ub = list(x = Inf))
norminc <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cinc),
  lb = list(x = 0), ub = list(x = Inf))
poispos <- glm(ypois ~ x, family = "poisson",
  method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = 0), ub = list(x = Inf))
poisinc <- glm(ypois ~ x, family = "poisson",
  method = cirls.fit, Cmat = list(x = cinc),
  lb = list(x = 0), ub = list(x = Inf))

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
  expect_equal(normpos$Cmat, cbind(0, cpos), ignore_attr = T)
  expect_equal(norminc$Cmat, cbind(0, cinc), ignore_attr = T)
  expect_equal(poispos$Cmat, cbind(0, cpos), ignore_attr = T)
  expect_equal(poisinc$Cmat, cbind(0, cinc), ignore_attr = T)

  # bounds are well computed
  expect_equal(normpos$lb, rep(0, p), ignore_attr = T)
  expect_equal(norminc$lb, rep(0, p - 1), ignore_attr = T)
  expect_equal(poispos$lb, rep(0, p), ignore_attr = T)
  expect_equal(poisinc$lb, rep(0, p - 1), ignore_attr = T)
  expect_equal(normpos$ub, rep(Inf, p), ignore_attr = T)
  expect_equal(norminc$ub, rep(Inf, p - 1), ignore_attr = T)
  expect_equal(poispos$ub, rep(Inf, p), ignore_attr = T)
  expect_equal(poisinc$ub, rep(Inf, p - 1), ignore_attr = T)

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
  lb = list(x = -Inf), ub = list(x = 0))
norminc_rev <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = -cinc),
  lb = list(x = -Inf), ub = list(x = 0))

# Check this is equivalent to previous constraint
test_that("reverting ub and lb is equal", {
  expect_mapequal(coef(normpos), coef(normpos_rev))
  expect_mapequal(coef(norminc), coef(norminc_rev))
})

# Equality constraint
betasum <- sum(betas)
normeq <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = list(x = betasum), ub = list(x = betasum))
test_that("equality constraint can be passed", {
  expect_equal(sum(coef(normeq)[-1]), betasum)
})

#----- Alternative solvers ------

# quadprog
quadprog_pos <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = 0), ub = list(x = Inf), qp_solver = "quadprog")
quadprog_inc <- glm(ypois ~ x, method = cirls.fit, Cmat = list(x = cinc),
  lb = list(x = 0), ub = list(x = Inf), qp_solver = "quadprog")
betasum <- sum(betas)
quadprog_eq <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = list(x = betasum), ub = list(x = betasum), qp_solver = "quadprog")

test_that("quadprog solver is integrated", {
  expect_true(all(coef(quadprog_pos)[-1] >= (0 - 1e-6)))
  expect_true(all(diff(coef(quadprog_inc)[-1]) >= (0 - 1e-6)))
  expect_equal(sum(coef(quadprog_eq)[-1]), betasum)
})

# coneproj
cone_pos <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = 0), ub = list(x = Inf), qp_solver = "coneproj")
cone_inc <- glm(ypois ~ x, method = cirls.fit, Cmat = list(x = cinc),
  lb = list(x = 0), ub = list(x = Inf), qp_solver = "coneproj")
betasum <- sum(betas)
cone_eq <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = list(x = betasum), ub = list(x = betasum), qp_solver = "coneproj")

test_that("coneproj solver is integrated", {
  expect_true(all(coef(cone_pos)[-1] >= (0 - 1e-6)))
  expect_true(all(diff(coef(cone_inc)[-1]) >= (0 - 1e-6)))
  expect_equal(sum(coef(cone_eq)[-1]), betasum)
})



#----- Test unconstrained model

# Apply base GLM
normglm <- glm(ynorm ~ x)
poisglm <- glm(ypois ~ x, family = "poisson")

# Apply CIRLS with no constraint (putting Inf)
normuncons <- glm(ynorm ~ x, method = cirls.fit, Cmat = list(x = cinc),
  lb = list(x = -Inf), ub = list(x = Inf))
poisuncons <- glm(ypois ~ x, family = "poisson",
  method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = -Inf), ub = list(x = Inf))

# Check they are identical
checkcomp <- c("coefficients", "residuals", "fitted.values", "effects", "R",
  "rank", "family", "linear.predictors", "deviance", "aic", "null.deviance",
  "weights", "df.residuals", "df.null")
test_that("Unconstrained model return the same as base GLM", {
  expect_equal(normglm[checkcomp], normuncons[checkcomp])
  expect_equal(poisglm[checkcomp], poisuncons[checkcomp])
})
