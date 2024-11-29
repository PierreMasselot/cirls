################################################################################
#
# Test that the cons_shape method works
#
################################################################################

# Libraries containing splines
library(splines)
library(dlnm)

#----- Generate data ------

# Parameters
n <- 1000
nsim <- 50

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
  Cmat <- cons_shape(basis, diff = 1)

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat)) |>
      predict())

  # Test results
  test_that(paste0("Monotone increasing constraints work with ", b), {
    expect_true(all(diff(res) >= -10^-10))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Increasing: ", b))
  lines(X, eta, lwd = 2)
}

#----- Monotone decreasing

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- cons_shape(basis, diff = 1, sign = -1)

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat)) |>
      predict())

  # Test results
  test_that(paste0("Monotone decreasing constraints work with ", b), {
    expect_true(all(diff(res) <= 10^-10))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Decreasing: ", b))
  lines(X, eta, lwd = 2)
}


#----- Convex

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- cons_shape(basis, diff = 2)

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat)) |>
      predict())

  # Test results
  test_that(paste0("Convex constraints work with ", b), {
    expect_true(all(diff(res, diff = 2) >= -10^-10))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Convex: ", b))
  lines(X, eta, lwd = 2)
}


#----- Concave

for (b in tested_bases){

  # Basis
  basis <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat <- cons_shape(basis, diff = 2, sign = -1)

  # Fit model
  res <- apply(Y, 2, function(y) glm(y ~ basis, family = "quasipoisson",
    method = "cirls.fit", Cmat = list(basis = Cmat)) |>
      predict())

  # Test results
  test_that(paste0("Concave constraints work with ", b), {
    expect_true(all(diff(res, diff = 2) <= 10^-10))
  })
  matplot(X, res, type = "l", col = "grey", main = paste0("Concave: ", b))
  lines(X, eta, lwd = 2)
}

#-------------------------
# Test the method for onebasis
#-------------------------

# Test that we get the right constraint matrix
test_that("We get the right constraint matrix with one basis ", {
for (b in tested_bases){

  # Basis
  basis1 <- onebasis(X, fun = b, df = p)
  basis2 <- do.call(b, list(x = X, df = p))

  # Generate constraint matrix
  Cmat1 <- cons_shape(basis1, diff = 1)
  Cmat2 <- cons_shape(basis2, diff = 1)

  # Test equality
  expect_true(all.equal(Cmat1, Cmat2))
}})

# Test for basis with no method
test_that("We get the default method for unknown onebasis", {
  strbasis <- onebasis(X, fun = "strata", df = p)
  expect_warning(Cmat <- cons_shape(strbasis))
  expect_true(all.equal(Cmat, diag(p)))
})
