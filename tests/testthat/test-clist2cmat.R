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

#----- Perform tests
test_that("it always returns a matrix", {
  expect_length(dim(clist2cmat(list(z = cz), mt, x)), 2)
})

test_that("matrix is expanded correctly",{
  # Only basis constrained
  expect_equal(
    clist2cmat(list(basis = cbas), mt, x),
    cbind(0, cbas, 0, 0, 0)
  )

  # Both basis and basis1
  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1), mt, x),
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0))
  )

  # Everything
  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1, z = cz), mt, x),
    rbind(cbind(0, cbas, 0, 0, 0), cbind(0, 0, 0, 0, 0, 0, cbas1, 0),
      c(rep(0, 8), 1))
  )
})

test_that("order of terms does not matter", {
  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1), mt, x),
    clist2cmat(list(basis1 = cbas1, basis = cbas), mt, x)
  )

  expect_equal(
    clist2cmat(list(basis = cbas, basis1 = cbas1, z = cz), mt, x),
    clist2cmat(list(basis1 = cbas1, z = cz, basis = cbas), mt, x)
  )
})

test_that("fails when list wringly specified", {
  # When names are wrongly specified
  expect_error(clist2cmat(list(x = cz), mt, x))
  expect_error(clist2cmat(list(basis12 = cbas1), mt, x))

  # When column number don't match
  expect_error(clist2cmat(list(basis = cbas1), mt, x))
  expect_error(clist2cmat(list(basis = cz), mt, x))
  expect_error(clist2cmat(list(basis1 = cbas), mt, x))
  expect_error(clist2cmat(list(basis1 = cz), mt, x))
  expect_error(clist2cmat(list(z = cbas), mt, x))
  expect_error(clist2cmat(list(z = cbas1), mt, x))
})
