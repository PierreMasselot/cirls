################################################################################
#
# Test that chack_cmat fills its role
#
################################################################################

#----- Underlying equality constraint

# Three parameters sum to zero
cmat1 <- rbind(rep(1, 3), rep(-1, 3))
check1 <- check_cmat(cmat1)

# Another constraint: both parameters non-negative and their sum nonpositive
# Indicates both are contrained to zero
cmat2 <- rbind(diag(2), -1)
check2 <- check_cmat(cmat2)

# Check it is detected
test_that("check_cmat detects equality constraints", {
  expect_equal(check1$equality, c(1, 2))
  expect_equal(check2$equality, 1:3)
})

#----- Reducible constraints

p <- 5

# Increasing convex
# because of convexity, the positive difference between coef 2 and 3 on is redundant
cmat3 <- rbind(diff(diag(p)), diff(diag(p), diff = 2))
check3 <- check_cmat(cmat3)

# S shape: although not of full row rank, it is irreducible
cmat4 <- rbind(
  diag(p)[1,], # positive
  diff(diag(p))[c(1, p - 1),], # Increasing at both end
  diff(diag(p), diff = 2)[1:(p/2 - 1),], # First half convex
  -diff(diag(p), diff = 2)[(p/2):(p-2),] # second half concave
)
check4 <- check_cmat(cmat4)

# Test
test_that("check_cmat takes good decisions on irreducibility", {
  expect_gt(length(check3$redundant), 0)
  expect_length(check4$redundant, 0)
})
