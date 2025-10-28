################################################################################
#
# Test functions when Cmat is not full row rank
#
################################################################################

#----- Parameters

# Sample size
n <- 500

# Number of predictors
p <- 10

#----- Generate

# X matrix
X <- matrix(rnorm(n * p), n, p)

# Coefficients
Beta <- exp((1:p) - (p / 2)) / (1 + exp((1:p) - (p / 2)))

# Response
Y <- X %*% Beta + rnorm(n)

#----- Fit model

# Constraint matrix: S-shape
Cmat <- rbind(
  diag(p)[1,], # positive
  diff(diag(p))[c(1, p - 1),], # Increasing at both end
  diff(diag(p), diff = 2)[1:(p/2 - 1),], # First half convex
  -diff(diag(p), diff = 2)[(p/2):(p-2),] # second half concave
)

# Fit model (no intercept to make sure Cmat has more lines than columns)
res <- glm(Y ~ 0 + X, method = "cirls.fit", Cmat = Cmat,
  lb = 0, ub = Inf)

# Check that model fits and is feasible
test_that("model fits well", {
  expect_false(any(is.na(coef(res))))
  expect_true(all((Cmat %*% coef(res)) >= (0 - 1e-6)))
})

#----- Check methods behave well
test_that("methods return warning",{
  expect_warning(vcov(res))
  expect_warning(confint(res))
})

v <- suppressWarnings(vcov(res))
ci <- suppressWarnings(confint(res))
test_that("methods return NAs", {
  expect_true(all(is.na(v)))
  expect_true(all(is.na(ci)))
})
