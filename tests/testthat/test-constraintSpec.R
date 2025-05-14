##################################################################
#
# test different possibilities to pass constraints
#
##################################################################

library(splines)

#----- Parameters
n <- 100

#----- Create objects

# Basis object
basis <- ns(rnorm(n), df = 5)
basis1 <- matrix(rnorm(n * 2), ncol = 2)

# Some covariate
z <- rnorm(n)
w <- rnorm(n)

# Some response
y <- rnorm(n)

# constraint matrices
cbas <- diff(diag(5)) # Increase
cbas2 <- diff(diag(5), diff = 2)
cbas1 <- diag(2) # All positive
cbasf <- shapeConstr(basis, "inc")
cz <- 1 # positive and not matrix object

# Terms and model matrix
lmres <- lm(y ~ basis + basis1 + z, x = T)
mf <- model.frame(lmres)

# control <- list(constr = ~ shape(basis1, "inc") + zerosum(basis), Cmat = list(z = cz))

#----------------------
# Test formula
#----------------------

#----- General formula

test_that("simple matrices work", {
  expect_identical(shapeConstr(basis, "inc"),
    constr2Clist(~ shape(basis, "inc"), mf = mf)[[1]],
    ignore_attr = "vars")
  expect_identical(shapeConstr(basis, "cvx"),
    constr2Clist(~ shape(basis, "cvx"), mf = mf)[[1]],
    ignore_attr = "vars")
  expect_identical(zerosumConstr(basis),
    constr2Clist(~ zerosum(basis), mf = mf)[[1]],
    ignore_attr = "vars")
})

test_that("two different formulations are identical", {
  expect_identical(constr2Clist(~ shapeConstr(basis, "inc"), mf = mf),
    constr2Clist(~ shape(basis, "inc"), mf = mf),
    ignore_attr = "names")
  expect_identical(constr2Clist(~ zerosumConstr(basis), mf = mf),
    constr2Clist(~ zerosum(basis), mf = mf),
    ignore_attr = "names")
})

test_that("it removes unknown function", {
  constr <- ~ shape(basis, "inc") + unknown(basis)
  expect_warning(res <- constr2Clist(constr, mf))
  expect_identical(res, constr2Clist(~ shape(basis, "inc"), mf))
})

test_that("it removes unknown terms", {
  constr <- ~ shape(basis, "inc") + zerosum(unknown)
  expect_warning(res <- constr2Clist(constr, mf))
  expect_identical(res, constr2Clist(~ shape(basis, "inc"), mf))
})

test_that("it works with more complex cases", {
  expect_identical(
    constr2Clist(~ shape(basis, "inc") + shape(basis, "cvx"), mf),
    list(shapeConstr(basis, "inc"), shapeConstr(basis, "cvx")),
    ignore_attr = c("names", "vars"))
  expect_identical(
    constr2Clist(~ shape(basis, "inc") + zerosum(basis), mf),
    list(shapeConstr(basis, "inc"), zerosumConstr(basis)),
    ignore_attr = c("names", "vars"))
})

form <- ~ shape(basis, "inc")
test_that("it works without the model frame", {
  expect_identical(shapeConstr(basis, "inc"),
    constr2Clist(form)[[1]],
    ignore_attr = "vars")
})


#----- Test transformation of list into matrix
test_that("it always returns a matrix", {
  expect_length(dim(buildCmat(Cmat = list(z = cz), mf = mf)$Cmat), 2)
})

test_that("matrix is expanded correctly",{
  # Only basis constrained
  expect_equal(
    buildCmat(Cmat = list(basis = cbas), mf = mf)$Cmat,
    cbind(0, cbas, 0, 0, 0),
    ignore_attr = T
  )

  # Both basis and basis1
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1), mf = mf)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0)),
    ignore_attr = T
  )

  # Everything
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1, z = cz), mf = mf)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0),
      c(rep(0, 8), 1)),
    ignore_attr = T
  )
})

test_that("fails when list wrongly specified", {
  # When names are wrongly specified
  expect_warning(buildCmat(Cmat = list(x = cz), mf = mf),
      regexp = "Unknown terms in Cmat") |>
    expect_warning(regexp = "No valid constraint")
  expect_warning(buildCmat(Cmat = list(basis12 = cbas1), mf = mf),
    regexp = "Unknown terms in Cmat") |>
    expect_warning(regexp = "No valid constraint")
  expect_warning(buildCmat(Cmat = list(basis1 = cbas1, x = cz), mf = mf),
    regexp = "Unknown terms in Cmat")

  # When column number don't match
  expect_error(buildCmat(Cmat = list(basis = cbas1), mf = mf),
    regexp = "inconsistent with model frame")
  expect_error(buildCmat(Cmat = list(basis = cz), mf = mf),
    regexp = "inconsistent with model frame")
  expect_error(buildCmat(Cmat = list(basis1 = cbas), mf = mf),
    regexp = "inconsistent with model frame")
  expect_error(buildCmat(Cmat = list(basis1 = cz), mf = mf),
    regexp = "inconsistent with model frame")
  expect_error(buildCmat(Cmat = list(z = cbas), mf = mf),
    regexp = "inconsistent with model frame")
  expect_error(buildCmat(Cmat = list(z = cbas1), mf = mf),
    regexp = "inconsistent with model frame")
})

test_that("buildCmat works with list instead of model.frame", {
  ml <- list(basis = basis, basis1 = basis1, z = z)

  expect_equal(
    buildCmat(Cmat = list(basis = cbas), mf = ml)$Cmat,
    cbind(0, cbas, 0, 0, 0),
    ignore_attr = T
  )

  # Both basis and basis1
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1), mf = ml)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0)),
    ignore_attr = T
  )

  # Everything
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1, z = cz), mf = ml)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0),
      c(rep(0, 8), 1)),
    ignore_attr = T
  )
})

#----- Check that not providing any Cmat results in same fit as glm

test_that("unconstrained model works", {

  # When nothing is provided
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = NULL, constr = NULL),
    "No valid constraint")
  expect_equal(coef(resu), coef(rescu))
  expect_equal(resu$deviance, rescu$deviance)

  # When all constraints are removed
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = list(x = cz)),
    "Unknown terms") |>
    expect_warning("No valid constraint")
  expect_equal(coef(resu), coef(rescu))
  expect_equal(resu$deviance, rescu$deviance)

  # With OSQP
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = NULL, constr = NULL, qp_solver = "osqp"),
    "No valid constraint")
  expect_equal(coef(resu), coef(rescu))
  expect_equal(resu$deviance, rescu$deviance)

  # With coneproj
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = NULL, constr = NULL, qp_solver = "coneproj"),
    "No valid constraint") |>
    expect_error("another solver")
})

#----- Test zero-sum constraint

# Build a zum to zero constraints
cm0 <- t(rep(1, ncol(basis) + ncol(basis1)))
attributes(cm0) <- c(attributes(cm0), list(lb = 0, ub = 0))
cm1 <- zerosumConstr(basis, basis1)
cm2 <- zerosumConstr(cbind(basis, basis1))

# Build a grouped one
cmg0 <- matrix(c(rep(1, ncol(basis)), rep(0, ncol(basis1)),
  rep(0, ncol(basis)), rep(1, ncol(basis1))), nrow = 2, byrow = TRUE)
attributes(cmg0) <- c(attributes(cmg0), list(lb = c(0, 0), ub = c(0, 0)))
cmg1 <- zerosumConstr(basis, basis1, group = TRUE)

test_that("zerosum works with severla terms", {
  expect_identical(cm0, cm1)
  expect_identical(cm0, cm2)
  expect_identical(cmg0, cmg1)
})

#----- Test complex constraint specification

# Build zero sum constraint with two terms
res1 <- buildCmat(mf = mf, constr = ~ zerosum(basis, basis1))
res2 <- buildCmat(mf = mf, Cmat = list("basis;basis1" = cm1),
  lb = list("basis;basis1" = attr(cm1, "lb")),
  ub = list("basis;basis1" = attr(cm1, "ub")))

test_that("Providing several terms work", {
  expect_identical(res1, res2, ignore_attr = TRUE)
  expect_identical(res1$Cmat, cbind(0, cm1, 0), ignore_attr = TRUE)
})

# Test bounds
test_that("Various ways to provide bounds work", {
  expect_identical(
    buildCmat(mf = mf, Cmat = list(basis = cbas)),
    buildCmat(mf = mf, Cmat = list(basis = cbas),
      lb = list(basis = 0), ub = list(basis = Inf))
  )
})
