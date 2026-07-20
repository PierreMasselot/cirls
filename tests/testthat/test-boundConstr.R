################################################################################
#
# Testing bound constraints
#
################################################################################

# Libraries containing splines
library(splines)
suppressMessages(library(dlnm))

#----- Generate data ------

# Parameters
n <- 1000
nsim <- 5

# Generate a sinusoid
X <- seq(0, 1, length.out = n)
eta <- .6 - X + (1.2 - X) * sin(X * 3 * pi)

# Generate several
set.seed(5)
Y <- rnorm(n, eta)

#-------------------------
# bs
#-------------------------

test_that("bounds work with bs", {

# Create basis
basis <- bs(X, df = 10)

# Different values (the default on the righ)
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE)

# Now on the left: this doesn't work
expect_warning(m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left")))

# This works
basis <- bs(X, df = 10, int = T)
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left"))
expect_equal(predict(m)[1], 1, ignore_attr = TRUE)

})

#-------------------------
# ps
#-------------------------

test_that("bounds work with ps", {

# Create basis
basis <- ps(X, df = 10)

# Different values (the default on the righ)
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

# sm cannot be reduced here
expect_error(m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, sm = 1)))

# Now on the left: this doesn't work
expect_warning(m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left")))

# This works
basis <- ps(X, df = 10, int = T)
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left"))
expect_equal(predict(m)[1], 1, ignore_attr = TRUE, tolerance = 10e-8)

})

#-------------------------
# ns
#-------------------------

test_that("bounds work with ns", {

# Create basis
basis <- ns(X, df = 10)

# Different values (the default on the righ)
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

# sm can be reduced here
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, sm = 1))
pred <- basis %*% coef(m)[-1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

# On the left
basis <- ns(X, df = 10, int = T)
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left"))
expect_equal(predict(m)[1], 1, ignore_attr = TRUE, tolerance = 10e-8)

# Both
basis <- ns(X, df = 10, int = T)
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 0, side = "both"))
pred <- predict(m)
expect_equal(pred[1], 0, ignore_attr = TRUE, tolerance = 10e-8)
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

})

#-------------------------
# factor
#-------------------------

test_that("bounds work for a factor", {

# Create basis
basis <- cut(X, breaks = 10)

# Different values (the default on the righ)
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

# Now on the left: this throws a warning and doesn't constrain
expect_warning(m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 1, side = "left")))

# This works
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left"))
expect_equal(predict(m)[1], 1, ignore_attr = TRUE, tolerance = 10e-8)

# Both sides
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "both"))
pred <- predict(m)
expect_equal(pred[1], 1, ignore_attr = TRUE, tolerance = 10e-8)
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

})

#-------------------------
# strata
#-------------------------

test_that("bounds work for a strata object", {

# Create basis
basis <- dlnm:::strata(X, df = 10)

# Different values (the default on the righ)
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

# Now on the left: this doesn't work
expect_warning(m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 1, side = "left")))

# This works
basis <- dlnm:::strata(X, df = 10, intercept = T)
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left"))
expect_equal(predict(m)[1], 1, ignore_attr = TRUE, tolerance = 10e-8)

# Both sides
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "both"))
pred <- predict(m)
expect_equal(pred[1], 1, ignore_attr = TRUE, tolerance = 10e-8)
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

})

#-------------------------
# onebasis
#-------------------------

test_that("bounds work for a onebasis object", {

# Create basis
basis <- onebasis(X, fun = "bs", df = 10)

# Different values (the default on the righ)
m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 0, sm = 1))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

m <- glm(Y ~ basis, method = "cirls.fit",
  constr = ~ bound(basis, value = 1))
pred <- predict(m) - coef(m)[1]
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

# Now on the left: this doesn't work
expect_warning(m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 1, side = "left")))

# This works
basis <- onebasis(X, fun = "bs", df = 10, intercept = TRUE)
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "left"))
expect_equal(predict(m)[1], 1, ignore_attr = TRUE, tolerance = 10e-8)

# Both sides
m <- glm(Y ~ basis - 1, method = "cirls.fit",
  constr = ~ bound(basis, value = 1, side = "both"))
pred <- predict(m)
expect_equal(pred[1], 1, ignore_attr = TRUE, tolerance = 10e-8)
expect_equal(tail(pred, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

})

#-------------------------
# Longer bounds
#-------------------------

test_that("Longer bounds work", {

  # B-splines
  basis <- bs(X, df = 10)
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 0, thr = .7))
  pred <- basis %*% coef(m)[-1]
  expect_all_equal(pred[X >= .7], 0)

  # Another value
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 1, thr = .7))
  pred <- basis %*% coef(m)[-1]
  expect_all_equal(pred[X >= .7], 1)

  # P-splines
  basis <- ps(X, df = 10)
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 0, thr = .7))
  pred <- basis %*% coef(m)[-1]
  expect_all_equal(pred[X >= .7], 0)

  # Natural splines
  basis <- ns(X, df = 10)
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 0, thr = .7))
  pred <- basis %*% coef(m)[-1]
  expect_all_equal(pred[X >= .7], 0)

  # Natural splines and another value
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 1, thr = .7))
  pred <- basis %*% coef(m)[-1]
  expect_all_equal(pred[X >= .7], 1)

  # Onebasis
  basis <- onebasis(X, fun = "bs", df = 10)
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 0, thr = .7))
  pred <- basis %*% coef(m)[-1]
  expect_all_equal(pred[X >= .7], 0)

  # Factor
  basis <- cut(X, breaks = 10)
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 0, thr = 8))
  pred <- predict(m) - coef(m)[1]
  expect_all_equal(pred[X >= .7], 0)

  # Another value
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 1, thr = 8))
  pred <- predict(m) - coef(m)[1]
  expect_all_equal(pred[X >= .7], 1)

  # Strata
  basis <- dlnm:::strata(X, df = 10)
  m <- glm(Y ~ basis, method = "cirls.fit",
    constr = ~ bound(basis, value = 0, thr = .7))
  pred <- predict(m) - coef(m)[1]
  expect_all_equal(pred[X >= .7], 0)
})


# plot(X, pred)
# abline(v = attr(basis, "knots"))
# abline(h = 0)
