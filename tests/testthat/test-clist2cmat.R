##################################################################
#
# test different possibilities to pass constraints
#
##################################################################

#----- Parameters
n <- 100

#----- Create objects

# Basis object
basis <- matrix(rnorm(n * 5), ncol = 5)
basis1 <- matrix(rnorm(n * 2), ncol = 2)

# Some covariate
z <- rnorm(n)

# Some response
y <- rnorm(n)

# constraint matrices
cbas <- diff(diag(5)) # Increase
cbas1 <- diag(2) # All positive
cz <- 1 # positive and not matrix object

# Terms and model matrix
lmres <- lm(y ~ basis + basis1 + z, x = T)
mt <- terms(lmres)
x <- lmres$x

#----- Test trasnformation of list into matrix
test_that("it always returns a matrix", {
  expect_length(dim(clist2cmat(list(z = cz), mt, x)), 2)
})

test_that("matrix is expanded correctly",{
  # Only basis constrained
  expect_equal(
    clist2cmat(list(basis = cbas), mt, x),
    cbind(0, cbas, 0, 0, 0),
    ignore_attr = T
  )

  # Both basis and basis1
  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1), mt, x),
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0)),
    ignore_attr = T
  )

  # Everything
  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1, z = cz), mt, x),
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0),
      c(rep(0, 8), 1)),
    ignore_attr = T
  )
})

test_that("order of terms does not matter", {
  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1), mt, x),
    clist2cmat(list(basis1 = cbas1, basis = cbas), mt, x),
    ignore_attr = T
  )

  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1, z = cz), mt, x),
    clist2cmat(list(basis1 = cbas1, z = cz, basis = cbas), mt, x),
    ignore_attr = T
  )
})

test_that("fails when list wrongly specified", {
  # When names are wrongly specified
  expect_error(expect_warning(clist2cmat(list(x = cz), mt, x)))
  expect_error(expect_warning(clist2cmat(list(basis12 = cbas1), mt, x)))
  expect_warning(clist2cmat(list(basis1 = cbas1, x = cz), mt, x))

  # When column number don't match
  expect_error(clist2cmat(list(basis = cbas1), mt, x))
  expect_error(clist2cmat(list(basis = cz), mt, x))
  expect_error(clist2cmat(list(basis1 = cbas), mt, x))
  expect_error(clist2cmat(list(basis1 = cz), mt, x))
  expect_error(clist2cmat(list(z = cbas), mt, x))
  expect_error(clist2cmat(list(z = cbas1), mt, x))
})

#----- Test with bounds

# Initialise list of constraints
Clist <- list(basis = cbas, basis1 = cbas1, z = cz)
Cmat <- clist2cmat(Clist, mt, x)

# Tests
test_that("works with bounds", {

  # Works as usual
  res <- glm(y ~ basis + basis1 + z, method = "cirls.fit", Cmat = Clist)
  expect_equal(Cmat, res$Cmat)

  # Works with vectors
  res <- glm(y ~ basis + basis1 + z, method = "cirls.fit", Cmat = Clist,
    lb = rep(0.5, 7))
  expect_equal(Cmat, res$Cmat)

  # Works with list
  res <- glm(y ~ basis + basis1 + z, method = "cirls.fit", Cmat = Clist,
    lb = list(z = 0.5, basis = 1))
  expect_equal(Cmat, res$Cmat)
  expect_equal(c(1, 1, 1, 1, 0, 0, .5), res$lb)

  # Issues warnings
  expect_warning(res2 <- glm(y ~ basis + basis1 + z, method = "cirls.fit",
    Cmat = Clist, lb = list(z = 0.5, basis = 1, wrong = -1)))
  expect_equal(Cmat, res2$Cmat)
  expect_equal(c(1, 1, 1, 1, 0, 0, .5), res2$lb)
})
