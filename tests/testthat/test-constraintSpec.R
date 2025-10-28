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
mt <- terms(mf)
assign <- attr(stats::model.matrix(mt, mf), "assign")
# control <- list(constr = ~ shape(basis1, "inc") + zerosum(basis), Cmat = list(z = cz))

#----------------------
# Test formula
#----------------------

#----- Test the underlying constr2Clist function

test_that("constr2Clist builds simple matrices", {
  expect_identical(shapeConstr(basis, "inc"),
    constr2Clist(str2lang("shape(basis, 'inc')"), label = "", int = 0, mf = mf))
  expect_identical(shapeConstr(basis, "cvx"),
    constr2Clist(str2lang("shape(basis, 'cvx')"), label = "", int = 0, mf = mf))
  expect_identical(zerosumConstr(basis),
    constr2Clist(str2lang("zerosum(basis)"), label = "", int = 0, mf = mf))
})

test_that("constr2Clist accepts two different formulations", {
  expect_identical(
    constr2Clist(str2lang("shape(basis, 'inc')"), label = "", int = 0, mf = mf),
    constr2Clist(str2lang("shapeConstr(basis, 'inc')"), label = "", int = 0,
      mf = mf))
  expect_identical(
    constr2Clist(str2lang("zerosum(basis)"), label = "", int = 0, mf = mf),
    constr2Clist(str2lang("zerosumConstr(basis)"), label = "", int = 0,
      mf = mf))
})

#----- Test formula part of buildCmat

test_that("buildCmat stops for unknown function", {
  constr <- ~ shape(basis, "inc") + unknown(z)
  expect_error(buildCmat(mf, constr = constr), regexp = "z")
})

test_that("buildCmat removes unknown terms", {
  constr <- ~ shape(basis, "inc") + zerosum(unknown)
  expect_warning(res <- buildCmat(mf, constr = constr), regexp = "unknown")
})

test_that("it works with more complex cases", {
  cm1 <- shapeConstr(basis, "inc")
  cm2 <- zerosumConstr(basis1, z)
  trueCmat <- list(Cmat = cbind(0, rbind(cm1$Cmat, 0),
    rbind(matrix(0, nrow = nrow(cm1$Cmat), ncol = ncol(cm2$Cmat)), cm2$Cmat)),
    lb = c(cm1$lb, cm2$lb), ub = c(cm1$ub, cm2$ub))
  expect_identical(
    buildCmat(mf, constr = ~ shape(basis, "inc") + zerosum(basis1, z)),
    trueCmat, ignore_attr = TRUE)
})


#----- Test transformation of list into matrix
test_that("it always returns a matrix", {
  expect_length(dim(buildCmat(Cmat = list(z = cz), lb = list(z = 0),
    ub = list(z = Inf), mf = mf)$Cmat), 2)
})

test_that("matrix is expanded correctly",{
  # Only basis constrained
  expect_equal(
    buildCmat(Cmat = list(basis = cbas), lb = list(basis = 0),
      ub = list(basis = Inf), mf = mf)$Cmat,
    cbind(0, cbas, 0, 0, 0),
    ignore_attr = TRUE
  )

  # Both basis and basis1
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1),
      lb = list(basis = 0, basis1 = 0),
      ub = list(basis = Inf, basis1 = Inf), mf = mf)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0)),
    ignore_attr = TRUE
  )

  # Everything
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1, z = cz),
      lb = list(basis = 0, basis1 = 0, z = 0),
      ub = list(basis = Inf, basis1 = Inf, z = Inf),
      mf = mf)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0),
      c(rep(0, 8), 1)),
    ignore_attr = TRUE
  )
})

test_that("buildCmat fails when list wrongly specified", {
  # When names are wrongly specified
  expect_warning(buildCmat(Cmat = list(x = cz), mf = mf),
      regexp = ": x") |> expect_warning("No valid constraint")
  expect_warning(buildCmat(Cmat = list(basis12 = cbas1), mf = mf),
    regexp = "basis12") |> expect_warning("No valid constraint")
  expect_warning(buildCmat(Cmat = list(basis1 = cbas1, x = cz), mf = mf),
    regexp = ": x") |> expect_warning(regexp = "No `lb` found") |>
    expect_warning(regexp = "No `ub` found")

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
    ignore_attr = TRUE
  ) |>
    expect_warning(regexp = "No `lb` found") |>
    expect_warning(regexp = "No `ub` found")

  # Both basis and basis1
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1),
      lb = list(basis = 0, basis1 = 0),
      ub = list(basis = Inf, basis1 = Inf),
      mf = ml)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0)),
    ignore_attr = TRUE
  )

  # Everything
  expect_equal(
    buildCmat(Cmat = list(basis = cbas, basis1 = cbas1, z = cz),
      lb = list(basis = 0, basis1 = 0, z = 0),
      ub = list(basis = Inf, basis1 = Inf, z = Inf),
      mf = ml)$Cmat,
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0),
      c(rep(0, 8), 1)),
    ignore_attr = TRUE
  )
})

#----- Check that not providing any Cmat results in same fit as glm

test_that("unconstrained model works", {

  # When nothing is provided
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(
    rescu <- glm(y ~ basis + basis1 + z, method = "cirls.fit",
      Cmat = NULL, constr = NULL),
    "No constraint")
  expect_equal(coef(resu), coef(rescu))
  expect_equal(resu$deviance, rescu$deviance)

  # When all constraints are removed
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = list(x = cz)),
    "dropped term.*: x") |>
    expect_warning("No valid constraint")
  expect_equal(coef(resu), coef(rescu))
  expect_equal(resu$deviance, rescu$deviance)

  # With OSQP
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = NULL, constr = NULL, qp_solver = "osqp"),
    "No constraint")
  expect_equal(coef(resu), coef(rescu))
  expect_equal(resu$deviance, rescu$deviance)

  # With coneproj
  resu <- glm(y ~ basis + basis1 + z)
  expect_warning(rescu <- glm(y ~ basis + basis1 + z,
    method = "cirls.fit", Cmat = NULL, constr = NULL, qp_solver = "coneproj"),
    "No constraint") |>
    expect_error("another solver")
})

#----- Test zero-sum constraint

# Build a zum to zero constraints
cm0 <- list(Cmat = t(rep(1, ncol(basis) + ncol(basis1))), lb = 0, ub = 0)
cm1 <- zerosumConstr(basis, basis1)
cm2 <- zerosumConstr(cbind(basis, basis1))

# Build a grouped one
cmg0 <- list(Cmat = matrix(c(rep(1, ncol(basis)), rep(0, ncol(basis1)),
    rep(0, ncol(basis)), rep(1, ncol(basis1))), nrow = 2, byrow = TRUE),
  lb = c(0, 0), ub = c(0, 0))
cmg1 <- zerosumConstr(basis, basis1, group = TRUE)

test_that("zerosum works with several terms", {
  expect_identical(cm0, cm1)
  expect_identical(cm0, cm2)
  expect_identical(cmg0, cmg1)
})

#----- Test complex constraint specification

# Build zero sum constraint with two terms
res1 <- buildCmat(mf = mf, constr = ~ zerosum(basis, basis1))
res2 <- buildCmat(mf = mf,
  Cmat = list("basis;basis1" = t(rep(1, ncol(basis) + ncol(basis1)))),
  ub = list("basis;basis1" = 0), lb = list("basis;basis1" = 0))

test_that("Providing several terms work", {
  expect_identical(res1, res2, ignore_attr = TRUE)
  expect_identical(res1$Cmat,
    cbind(0, t(rep(1, ncol(basis) + ncol(basis1))), 0),
    ignore_attr = TRUE)
})


