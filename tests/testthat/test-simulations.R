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

# Test that we have the same as an unconstrained model called directly
checkcomp <- c("coefficients", "residuals", "fitted.values", "effects", "R",
  "rank", "family", "linear.predictors", "deviance", "aic", "null.deviance",
  "weights", "df.residuals", "df.null")
test_that("Unconstrained model return the same as base GLM", {
  expect_equal(res0[checkcomp], uc1[checkcomp])
  expect_equal(res0[checkcomp], uc2[checkcomp])
})

#----- Test with data within a data.frame
data <- data.frame(Y = y, X = x)
form <- sprintf("Y ~ %s", paste("X.", 1:p, sep = "", collapse = " + "))

# Unconstrained
udf <- glm(form, data = data)

# Predictor literal in formula
resdf <- glm(form, data = data, method = cirls.fit, Cmat = cbind(0, cinc))

# Predictor implied
resdot <- glm(Y ~ ., data = data, method = cirls.fit, Cmat = cbind(0, cinc))

test_that("`uncons` works with `data` argument", {
  expect_equal(udf[checkcomp], uncons(resdf)[checkcomp])
  expect_equal(udf[checkcomp], uncons(resdot)[checkcomp])
})

#------------------------------
# Simulate
#------------------------------

#----- With constraints

# Simulate
simcons <- simulCoef(res1, nsim = 1000, seed = 5)

# Plot
plot(simcons[,3:4], pch = ".", cex = 2)
abline(a = 0, b = 1, col = 2)

#----- Unconstrained

# Simulate
simun <- simulCoef(res1, nsim = 1000, seed = 5, cons = F)

# Plot
plot(simun[,3:4], pch = ".", cex = 2)
abline(a = 0, b = 1, col = 2)


