################################################################################
#
# Test that the shapeConstr method works
#
################################################################################

# Libraries containing splines
library(splines)
suppressMessages(library(dlnm))

#----- Generate data ------

# Parameters
n <- 1000
nsim <- 5

# Generate a simple parabolic relationship
#set.seed(1)
X <- seq(0, 1, length.out = n)
# eta <- log(10) + (x-.3)^2
eta <- log(10) - .5 * sin(X * 1.6 * pi)

# Generate several
set.seed(5)
Y <- replicate(nsim, rpois(n, exp(eta)))

#-------------------------
# Tests
#-------------------------

# Parameters
tested_bases <- c("ps", "bs", "ns")
p <- 10

#----- Monotone increasing

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- shapeConstr(basis, shape = "inc")

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat$Cmat),
    lb = list(basis = Cmat$lb), ub = list(basis = Cmat$ub)) |>
      predict())

  # Test results
  test_that(paste0("Monotone increasing constraints work with ", b), {
    expect_true(all(diff(res) >= -sqrt(.Machine$double.eps)))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Increasing: ", b))
  lines(X, eta, lwd = 2)
}

#----- Monotone decreasing

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- shapeConstr(basis, shape = "dec")

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat$Cmat),
    lb = list(basis = Cmat$lb), ub = list(basis = Cmat$ub)) |>
      predict())

  # Test results
  test_that(paste0("Monotone decreasing constraints work with ", b), {
    expect_true(all(diff(res) <= sqrt(.Machine$double.eps)))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Decreasing: ", b))
  lines(X, eta, lwd = 2)
}


#----- Convex

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- shapeConstr(basis, shape = "cvx")

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat$Cmat),
    lb = list(basis = Cmat$lb), ub = list(basis = Cmat$ub)) |>
      predict())

  # Test results
  test_that(paste0("Convex constraints work with ", b), {
    expect_true(all(diff(res, diff = 2) >= -sqrt(.Machine$double.eps)))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Convex: ", b))
  lines(X, eta, lwd = 2)
}


#----- Concave

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- shapeConstr(basis, shape = "ccv")

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat$Cmat),
    lb = list(basis = Cmat$lb), ub = list(basis = Cmat$ub)) |>
      predict())

  # Test results
  test_that(paste0("Concave constraints work with ", b), {
    expect_true(all(diff(res, diff = 2) <= sqrt(.Machine$double.eps)))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Concave: ", b))
  lines(X, eta, lwd = 2)
}

#-------------------------
# Test the method for onebasis
#-------------------------

# Test that we get the right constraint matrix
test_that("We get the right constraint matrix with `onebasis`", {
for (b in tested_bases){

  # Basis
  basis1 <- onebasis(X, fun = b, df = p)
  basis2 <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat1 <- shapeConstr(basis1, shape = "inc")
  Cmat2 <- shapeConstr(basis2, shape = "inc")

  # Test equality
  expect_true(all.equal(Cmat1, Cmat2))
}})

# Test for basis with no method
test_that("We get the default method for unknown `onebasis`", {
  strbasis <- onebasis(X, fun = "strata", df = p)
  expect_warning(Cmat <- shapeConstr(strbasis, shape = "pos"))
  expect_identical(Cmat$Cmat, diag(p), ignore_attr = TRUE)
})


#-------------------------
# Antonio's test
#-------------------------

# Generate data
x1 <- 1:100
mu <- x1/100 - (x1/100)^2 - 0.3*(x1/100)^3
set.seed(13041975)
y <- mu + rnorm(100,0,0.1)

# Spline bases
ks <- 1:4*20
X <- ns(x1, knots = ks)

# Fit
Cmat <- shapeConstr(X, shape = c("dec", "ccv"))
m <- glm(y ~ X, method = cirls.fit, Cmat = list(X=Cmat$Cmat),
  lb = list(X = Cmat$lb), ub = list(X = Cmat$ub))
pred <- predict(m)

# Formal test
test_that("More challengin test for ns", {
  expect_true(all(diff(pred, diff = 2) <= sqrt(.Machine$double.eps)))
  expect_true(all(diff(pred) <= sqrt(.Machine$double.eps)))
})

# Plot
plot(x1,y)
lines(x1,mu,type="l",col=1,lwd=2)
lines(x1,pred,lwd=2, col=3)

#-------------------------
# Factor
#-------------------------

# Parameters
n <- 1000
nsim <- 50

# Generate a simple parabolic relationship
X <- sample(1:10, n, replace = T)
# eta <- log(10) + (x-.3)^2
eta <- - .5 * sin(X * 1.6 * pi / 10)

# Generate several
Y <- rpois(n, exp(eta))

#----- Factor

# Transform as factor
Xf <- factor(X)

# Constraint matrix
Cmat <- shapeConstr(Xf, "inc")

# Basic test
treat <- glm(Y ~ Xf, family = "quasipoisson", method = "cirls.fit",
  Cmat = list(Xf = Cmat$Cmat), lb = list(Xf = Cmat$lb), ub = list(Xf = Cmat$ub))
plot(X, eta, pch = 16, ylim = range(c(eta, predict(treat))))
points(X, predict(treat), pch = 15, col = 3)

# When intercept is "included" in the factor
Cmat0 <- shapeConstr(Xf, "inc", intercept = TRUE)
int <- glm(Y ~ 0 + Xf, family = "quasipoisson", method = "cirls.fit",
  Cmat = list(Xf = Cmat0$Cmat), lb = list(Xf = Cmat0$lb),
  ub = list(Xf = Cmat0$ub))
plot(X, eta, pch = 16, ylim = range(c(eta, predict(int))))
points(X, predict(int), pch = 15, col = 3)

# Automatic addition of intercept
int2 <- glm(Y ~ 0 + Xf, family = "quasipoisson", method = "cirls.fit",
 constr = ~ shape(Xf, "inc"))
plot(X, eta, pch = 16, ylim = range(c(eta, predict(int2))))
points(X, predict(int2), pch = 15, col = 3)

# Helmert contrasts
Xf2 <- Xf
contrasts(Xf2) <- "contr.helmert"
Cmat2 <- shapeConstr(Xf2, "inc")
helm <- glm(Y ~ Xf2, family = "quasipoisson", method = "cirls.fit",
  Cmat = list(Xf2 = Cmat2$Cmat), lb = list(Xf2 = Cmat2$lb),
  ub = list(Xf2 = Cmat2$ub))
plot(X, eta, pch = 16, ylim = range(c(eta, predict(helm))))
points(X, predict(helm), pch = 15, col = 3)
# model.matrix(helm)

# Tests
test_that("Shape constraints on factors work", {
  expect_true(all(
    diff(predict(treat, newdata = data.frame(Xf = levels(Xf)))) >=
      -sqrt(.Machine$double.eps)))
  expect_true(all(
    diff(predict(int, newdata = data.frame(Xf = levels(Xf)))) >=
      -sqrt(.Machine$double.eps)))
  expect_true(all(
    diff(predict(int2, newdata = data.frame(Xf = levels(Xf)))) >=
      -sqrt(.Machine$double.eps)))
  expect_true(all(
    diff(predict(helm, newdata = data.frame(Xf2 = levels(Xf2)))) >=
      -sqrt(.Machine$double.eps)))
})


###############################################################################
# Constraints on subdomains
###############################################################################

# #---------------------
# # Test constraining specific domain
# #---------------------
#
# #----- Generate data
#
# n <- 1000
#
# X <- seq(0, 1, length.out = n)
# # eta <- log(10) + (x-.3)^2
# eta <- 5 - exp(X)
#
# # Generate several
# set.seed(6)
# Y <- eta + rnorm(n, sd = .5)
#
# # A P-spline basis
# basis <- ps(x = X, df = 10)
#
# # A sequence of breaks
# brseq <- seq(-0.1, 1.1, by = .05)
#
# test_that("subdomain constraints work", {
#
#   #----- Initial test with a non-decreasing constraint
#   for (br in brseq){
#
#     # Fit a non-decreasing curve above br
#     Cmat1 <- mkDmat(d = 1, s = 1, knots = attr(basis, "knots"),
#       ord = attr(x, "degree") + 1, intercept = attr(basis, "intercept"),
#       lower = br, upper = Inf)
#     curve1 <- glm(Y ~ basis, method = "cirls.fit", Cmat = list(basis = Cmat1)) |>
#       predict()
#
#     # Test it works
#     expect_true(all(!isFALSE(
#       diff(curve1[X > br]) >= -sqrt(.Machine$double.eps))))
#     expect_true(any(!isFALSE(
#       diff(curve1[X < br]) < -sqrt(.Machine$double.eps))))
#
#     # Fit a non-decreasing curve below br
#     Cmat2 <- mkDmat(d = 1, s = 1, knots = attr(basis, "knots"),
#       ord = attr(x, "degree") + 1, intercept = attr(basis, "intercept"),
#       lower = -Inf, upper = br)
#     curve2 <- glm(Y ~ basis, method = "cirls.fit", Cmat = list(basis = Cmat2)) |>
#       predict()
#
#     # Test it works
#     expect_true(all(!isFALSE(
#       diff(curve2[X < br]) >= -sqrt(.Machine$double.eps))))
#     expect_true(any(!isFALSE(
#       diff(curve2[X > br]) < -sqrt(.Machine$double.eps))))
#
#     # Plot the curves
#     plot(X, eta, type = "l")
#     lines(X, curve1, col = 3)
#     lines(X, curve2, col = 4)
#     abline(v = attr(basis, "knots"), lty = 3)
#     abline(v = br, col = 2)
#   }
#
#   #----- Initial test with convexity constraint
#
#   # Testing with first degree constraint
#   for (br in brseq){
#
#     # Fit convex curve below br
#     Cmat1 <- mkDmat(d = 2, s = 1, knots = attr(basis, "knots"),
#       ord = attr(x, "degree") + 1, intercept = attr(basis, "intercept"),
#       lower = br, upper = Inf)
#     curve1 <- glm(Y ~ basis, method = "cirls.fit", Cmat = list(basis = Cmat1)) |>
#       predict()
#
#     # Test it works
#     expect_true(all(!isFALSE(
#       diff(curve1[X > br], diff = 2) >= -sqrt(.Machine$double.eps))))
#     expect_true(any(!isFALSE(
#       diff(curve1[X < br], diff = 2) < -sqrt(.Machine$double.eps))))
#
#     # Fit convex curve above br
#     Cmat2 <- mkDmat(d = 2, s = 1, knots = attr(basis, "knots"),
#       ord = attr(x, "degree") + 1, intercept = attr(basis, "intercept"),
#       lower = -Inf, upper = br)
#     curve2 <- glm(Y ~ basis, method = "cirls.fit", Cmat = list(basis = Cmat2)) |>
#       predict()
#
#     # Test it works
#     expect_true(all(!isFALSE(
#       diff(curve2[X < br], diff = 2) >= -sqrt(.Machine$double.eps))))
#     expect_true(any(!isFALSE(
#       diff(curve2[X > br], diff = 2) < -sqrt(.Machine$double.eps))))
#
#     # Plot curves
#     plot(X, eta, type = "l")
#     lines(X, curve1, col = 3)
#     lines(X, curve2, col = 4)
#     abline(v = attr(basis, "knots"), lty = 3)
#     abline(v = br, col = 2)
#   }
#
# })
#
