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

#------------------------------
# Test uncons function
#------------------------------


#----- Fit models

# Unconstrained
res0 <- glm(y ~ x)

# Increasing
cinc <- diff(diag(p))

# Cmat as an argument
res1 <- glm(y ~ x, method = cirls.fit, Cmat = list(x = cinc))

# Cmat in control list
res2 <- glm(y ~ x, method = cirls.fit, control = list(Cmat = list(x = cinc)))

#----- Basic use of uncons

uc1 <- uncons(res1)
uc2 <- uncons(res2)

# Test they are the same as an unconstrained model called directly
test_that("`uncons` fits unconstrained model", {
  expect_identical(res0, uc1)
  expect_identical(res0, uc2)
})

#----- Within function

f <- function(object) coef(uncons(object))
f1 <- f(res1)
f2 <- f(res2)

# test
test_that("`uncons` can be used within another function", {
  expect_identical(coef(res0), f1)
  expect_identical(coef(res0), f2)
})
