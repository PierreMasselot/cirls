################################################################################
#
# Solver comparison
#
################################################################################

library(cirls)

#--------------------
# Simple monotone regression (from section 4.1 of Meyer 2013)
#--------------------

# Data generation
n <- 50
x <- seq(0, 1, length.out = n)
eta <- 1 - (1 - x)^2
y <- eta + rnorm(n)

# Constraint matrix
Cmat <- cbind(0, 1, c(0, 2))

# Fit models
quadprogres <- glm(y ~ x + I(x^2), method = "cirls.fit", Cmat = Cmat,
  qp_solver = "quadprog") |> predict()
osqpres <- glm(y ~ x + I(x^2), method = "cirls.fit", Cmat = Cmat,
  qp_solver = "osqp") |> predict()
coneprojres <- glm(y ~ x + I(x^2), method = "cirls.fit", Cmat = Cmat,
  qp_solver = "coneproj") |> predict()

# Plot
plot(x, y)
matlines(x, cbind(eta, quadprogres, osqpres, coneprojres))

#--------------------
# More challenging
#--------------------

set.seed(1)
# Data generation
n <- 1000
p <- 10
x <- replicate(p, runif(n))
beta <- runif(p, 0, .1) |> sort()
eta <- x %*% beta
y <- eta + rnorm(n, sd = 10)

# Constraint matrix
Cmat <- diff(diag(p + 1))

# Fit models
quadprogres <- glm(y ~ x, method = "cirls.fit", Cmat = Cmat,
  qp_solver = "quadprog") |> coef()
osqpres <- glm(y ~ x, method = "cirls.fit", Cmat = Cmat,
  qp_solver = "osqp") |> coef()
coneprojres <- glm(y ~ x, method = "cirls.fit", Cmat = Cmat,
  qp_solver = "coneproj") |> coef()

# Compare
plot(beta, pch = 16, xlab = "")
matpoints(cbind(quadprogres, osqpres, coneprojres)[-1,])
