#------------------------------
# Generate some data
#------------------------------

#----- Parameters

# Number of obs
n <- 1000

# Coefficients
betas <- c(1, 1)
p <- length(betas)

#----- Generate data

# Uniform values between 0 and 1
x <- matrix(rnorm(n * p), n, p)

# Linear predictor
eta <- x %*% betas

# Simulate responses
y <- eta + rnorm(n, 0, .2)

#----- Fit model
Cmat <- diff(diag(p))
res <- glm(y ~ 0 + x, method = cirls.fit, Cmat = list(x = Cmat))

#------------------------------
# Test uncons
#------------------------------

# Fit the same model but without any constraint
resu <- glm(y ~ 0 + x)

# Another test
w <- sample(1:5, n, replace = TRUE)
res1 <- glm(y ~ x, weights = w, method = cirls.fit,
  Cmat = list(x = diff(diag(p))))
res2 <- glm(y ~ x, weights = w)

# Compare
test_that("`uncons` works", {
  expect_identical(uncons(res), resu)
  expect_identical(uncons(res1), res2)
})

#----- Test with data within a data.frame
data <- data.frame(Y = y, X = x)

# Predictor literal in formula
resdf <- glm(Y ~ 0 + X.1 + X.2, data = data, method = cirls.fit, Cmat = Cmat)
resdfu <- glm(Y ~ 0 + X.1 + X.2, data = data)

# Predictor implied
resdot <- glm(Y ~ ., data = data, method = cirls.fit, Cmat = cbind(0, Cmat))
resdotu <- glm(Y ~ ., data = data)

test_that("`uncons` works with `data` argument", {
  expect_identical(uncons(resdf), resdfu)
  expect_identical(uncons(resdot), resdotu)
})

#------------------------------
# Simulate
#------------------------------

#----- With constraints

# Simulate
simcons <- simulCoef(res, nsim = 1000, seed = 5)

# Plot
plot(simcons)
abline(a = 0, b = 1, col = 2)

#----- Unconstrained

# Simulate
simun <- simulCoef(res, nsim = 1000, seed = 5, cons = F)

# Plot
plot(simun)
abline(a = 0, b = 1, col = 2)


