################################################################################
#
#  Simple tests
#
################################################################################

#------------------------------
# Gaussian test
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
poispos <- glm(ypois ~ x, method = cirls.fit, Cmat = list(x = cpos))
poisinc <- glm(ypois ~ x, method = cirls.fit, Cmat = list(x = cinc))

#----- Tests ------

test_that("output contains all new elements", {

  # Contains all new outputs
  expect_true(all(c("coef.unconstrained", "active.cons", "Cmat", "bvec") %in%
    names(normpos)))
  expect_true(all(c("coef.unconstrained", "active.cons", "Cmat", "bvec") %in%
    names(norminc)))
  expect_true(all(c("coef.unconstrained", "active.cons", "Cmat", "bvec") %in%
    names(poispos)))
  expect_true(all(c("coef.unconstrained", "active.cons", "Cmat", "bvec") %in%
    names(poisinc)))

  # Cmat is correct
  expect_equal(normpos$Cmat, cbind(0, cpos))
  expect_equal(norminc$Cmat, cbind(0, cinc))
  expect_equal(poispos$Cmat, cbind(0, cpos))
  expect_equal(poisinc$Cmat, cbind(0, cinc))

  # bvec is correct
  expect_equal(normpos$bvec, rep(0 + formals(cirls.control)$bvectol, p))
  expect_equal(norminc$bvec, rep(0 + formals(cirls.control)$bvectol, p - 1))
  expect_equal(poispos$bvec, rep(0 + formals(cirls.control)$bvectol, p))
  expect_equal(poisinc$bvec, rep(0 + formals(cirls.control)$bvectol, p - 1))

})

test_that("results respect constraints", {
  # positive
  expect_true(all(coef(normpos)[-1] > 0))
  expect_true(all(coef(poispos)[-1] > 0))

  # increases
  expect_true(all(diff(coef(norminc)[-1]) > 0))
  expect_true(all(diff(coef(poisinc)[-1]) > 0))
})
