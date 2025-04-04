################################################################################
#
#  Test main routine
#
################################################################################

#------------------------------
# Generate some data
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
y <- eta + rnorm(n, 0, .2)

#----- Fit model

# Increasing
cinc <- diff(diag(p))

# Fit model
res <- glm(y ~ x, method = cirls.fit, Cmat = list(x = cinc))

#------------------------------
# Test methods
#------------------------------

test_that("S3 methods work", {
  expect_true(identical(vcov(res, seed = 1), vcov.cirls(res, seed = 1)))
  expect_true(identical(confint(res, seed = 1), confint.cirls(res, seed = 1)))
})

#------------------------------
# Perform inference
#------------------------------

#----- Simulations

# Simulate coefficients
simcoef <- simulCoef(res, nsim = 100)

# Test they all respect constraint
test_that("simulated coefficients respect constraints", {
  cons <- cinc %*% t(simcoef[,-1])
  expect_true(all(cons >= 0))
})

#----- Variance covariance

# Get variance covariance
cvcov <- vcov.cirls(res)

# Check variance are all positive
test_that("vcov matrix is well defined", {
  expect_true(all(diag(cvcov) > 0))
})

#----- Confidence intervals

# Get CIs
set.seed(1); cci95 <- confint.cirls(res, level = .95)
set.seed(1); cci90 <- confint.cirls(res, level = .90)
set.seed(1); cci80 <- confint.cirls(res, level = .80)

# Check confidence intervals include estimated coefficients
test_that("CI include estimated", {
  expect_true(all(cci95[,1] <= coef(res) & cci95[,2] >= coef(res)))
  expect_true(all(cci90[,1] <= coef(res) & cci90[,2] >= coef(res)))
  expect_true(all(cci80[,1] <= coef(res) & cci80[,2] >= coef(res)))
})

# Check confidence intervals are of increasing length
test_that("CI are nested", {
  expect_true(all(cci95[,1] <= cci90[,1]))
  expect_true(all(cci90[,1] <= cci80[,1]))
  expect_true(all(cci95[,2] >= cci90[,2]))
  expect_true(all(cci90[,2] >= cci80[,2]))
})


#------------------------------
# Tests with equality constraints
#------------------------------

#----- Sum equal to specific number

# Fit model
betasum <- sum(betas)
reseq1 <- glm(y ~ x, method = cirls.fit, Cmat = list(x = t(rep(1, p))),
  lb = betasum, ub = betasum)

# Compute vcov and ci
v <- vcov(reseq1)
ci <- confint(reseq1)

# Test
test_that("Inference with sum equality works", {
  expect_equal(sum(is.na(v)), 0)
  expect_equal(dim(v), rep(p + 1, 2))
  expect_true(all(diag(v) >= 0))
  expect_equal(sum(is.na(ci)), 0)
  expect_equal(dim(ci), c(p + 1, 2))
  expect_true(all(ci[,1] <= ci[,2]))
})

#----- Constrain specific coefficient

# Fit model
reseq2 <- glm(y ~ x, method = cirls.fit,
  Cmat = list(x = t(c(1, rep(0, p - 1)))), ub = 0)

# Compute vcov and ci
v <- vcov(reseq2)
ci <- confint(reseq2)

# Test
test_that("Inference with equality constraint on coefficient works", {
  expect_equal(v[2,2], 0)
  expect_equal(ci[2,1], ci[2,1])
})
