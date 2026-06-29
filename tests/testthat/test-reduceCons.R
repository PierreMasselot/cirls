################################################################################
#
# Test that chack_cmat fills its role
#
################################################################################

library(splines)

#----- Redundant constraints

p <- 10

test_that("it reduces simple matrices", {

  # Increasing convex
  # because of convexity, the positive difference between coef 2 and 3 on is redundant
  cmat1 <- rbind(diff(diag(p)), diff(diag(p), diff = 2))
  expect_warning(check1 <- reduceCons(cmat1, warn = TRUE))
  expect_lt(NROW(check1$Cmat), NROW(cmat1))

  # S shape: although not of full row rank, it is irreducible
  cmat2 <- rbind(
    diag(p)[1,], # positive
    diff(diag(p))[c(1, p - 1),], # Increasing at both end
    diff(diag(p), diff = 2)[1:(p/2 - 1),], # First half convex
    -diff(diag(p), diff = 2)[(p/2):(p-2),] # second half concave
  )
  expect_no_warning(check2 <- reduceCons(cmat2, warn = TRUE))
  expect_equal(length(check2$redundant), 0)
})


#----- Some cases with actual equality constraints

test_that("it works in some specficic cases", {

  # A simple test in which one ineq is redundant because of eq
  Cmat <- t(matrix(rep(0:1, 2), 2, 2))
  lb <- rep(0, 2)
  ub <- c(0, Inf)
  chk <- reduceCons(Cmat, lb, ub)
  chk2 <- reduceCons(Cmat, lb, rev(ub))
  expect_identical(chk[c("Cmat", "lb", "ub")], chk2[c("Cmat", "lb", "ub")])

  # Another simple test
  Cmat <- t(matrix(rep(0:1, 2), 2, 2))
  lb <- c(1, 0)
  ub <- c(1, Inf)
  expect_warning(chk <- reduceCons(Cmat, lb, ub, warn = TRUE))
  expect_equal(chk$lb, 1)
  expect_equal(chk$ub, 1)

  # Two redudant constraints, including an equality one
  cmlist <- Map(function(x, y) rbind(as.matrix(x), as.matrix(y)),
    boundConstr(diag(5)), shapeConstr(diag(5), shape = "pos"))
  cmlist$warn <- TRUE
  expect_warning(chk <- do.call(reduceCons, cmlist))

})


test_that("it reduces weird ns constraints", {

  # ns basis
  x <- ns(0:10, df = 5)

  # shape and bound constraints
  # The two don't give the exact same values but agree on the number
  cm1 <- Map(function(x, y) rbind(as.matrix(x), as.matrix(y)),
    boundConstr(x), shapeConstr(x, shape = "pos"))
  cm1$warn <- TRUE
  expect_warning(do.call(reduceCons, cm1))


  cm2 <- Map(function(x, y) rbind(as.matrix(x), as.matrix(y)),
    shapeConstr(x, shape = "pos"), boundConstr(x))
  cm2$warn <- TRUE
  expect_warning(do.call(reduceCons, cm2))

})

#----- Underlying equality constraint

test_that("it reduces underlying equality constraints", {

  # Example of undelrying equality to zero for both variables
  Cmat <- rbind(c(1, 0), c(0, -1), diff(diag(2)))
  lb <- rep(0, 3)
  ub <- rep(Inf, 3)
  expect_warning(chk <- reduceCons(Cmat, lb, ub, warn = TRUE))
  expect_equal(abs(chk$Cmat), diag(2))
  expect_equal(chk$lb, rep(0, 2))
  expect_equal(chk$ub, rep(0, 2))

  # Slightly more complex, both equal to 1
  Cmat <- rbind(diag(2), diff(diag(2)))
  lb <- c(1, -Inf, 0)
  ub <- c(Inf, 1, Inf)
  expect_warning(chk <- reduceCons(Cmat, lb, ub, warn = TRUE))
  expect_equal(abs(chk$Cmat), diag(2))
  expect_equal(chk$lb, rep(1, 2))
  expect_equal(chk$ub, rep(1, 2))

  # Equality constraint written as two inequalities
  Cmat <- t(matrix(0:1, 2, 2))
  lb <- c(0, -Inf)
  ub <- c(Inf, 0)
  expect_warning(chk <- reduceCons(Cmat, lb, ub, warn = TRUE))
  expect_equal(drop(chk$Cmat), c(0,1))
  expect_equal(chk$lb, 0)
  expect_equal(chk$ub, 0)
})

#----- Test no reducing

test_that("the function doesn't break when everything switched off",{

  # Some redundant constraints
  cmat1 <- rbind(diff(diag(p)), diff(diag(p), diff = 2))

  # Check switching equality off
  noeq <- reduceCons(cmat1, equality = FALSE)
  expect_equal(NROW(noeq$Cmat), p - 1)

  # Check switching both off
  nored <- reduceCons(cmat1, equality = FALSE, redundant = FALSE)
  expect_equal(NROW(nored$Cmat), NROW(cmat1))
})
