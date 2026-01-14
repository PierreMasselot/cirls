################################################################################
#
# Test that chack_cmat fills its role
#
################################################################################

#----- Underlying equality constraint

test_that("checkCmat detects equality constraints", {

  # Three parameters sum to zero
  cmat1 <- rbind(rep(1, 3), rep(-1, 3))
  expect_warning(check1 <- checkCmat(cmat1))
  expect_equal(check1$equality, rep(1, 2))

  # Another constraint: both parameters non-negative and their sum nonpositive
  # Indicates both are constrained to zero
  cmat2 <- rbind(diag(2), -1)
  expect_warning(check2 <- checkCmat(cmat2))
  expect_equal(check2$equality, rep(1, 3))

  # Two different equality constraints
  cmat3 <- rbind(cbind(cmat1, matrix(0, 2, 2)),
    cbind(matrix(0, 3, 3), cmat2))
  expect_warning(check3 <- checkCmat(cmat3))
  expect_equal(check3$equality, rep(1:2, c(2, 3)))
})



#----- Reducible constraints

p <- 10

test_that("checkCmat takes good decisions on irreducibility", {

  # Increasing convex
  # because of convexity, the positive difference between coef 2 and 3 on is redundant
  cmat1 <- rbind(diff(diag(p)), diff(diag(p), diff = 2))
  expect_warning(check1 <- checkCmat(cmat1, reduce = TRUE))
  expect_gt(sum(check1$redundant), 0)
  expect_lt(NROW(check1$Cmat), NROW(cmat1))

  # Removing the reduction
  expect_equal(cmat1, checkCmat(cmat1, reduce = FALSE, warn = FALSE)$Cmat)

  # S shape: although not of full row rank, it is irreducible
  cmat2 <- rbind(
    diag(p)[1,], # positive
    diff(diag(p))[c(1, p - 1),], # Increasing at both end
    diff(diag(p), diff = 2)[1:(p/2 - 1),], # First half convex
    -diff(diag(p), diff = 2)[(p/2):(p-2),] # second half concave
  )
  expect_no_warning(check2 <- checkCmat(cmat2))
  expect_equal(sum(check2$redundant), 0)
})
