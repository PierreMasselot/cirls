################################################################################
#
# Testing that the main function runs well
#
################################################################################

#----- Generate data

# Parameters
n <- 1000
betas <- c(0, 1, 2, -1, 1)
p <- length(betas)

# Generate
set.seed(5678)
x <- matrix(rnorm(n * p), n, p)
eta <- 1 + x %*% betas
ynorm <- eta + rnorm(n, 0, .2)

#-----------------------
# Try to call with do.call
#-----------------------

# argument list
cpos <- diag(p)
arglist <- list(formula = ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = 0), ub = list(x = Inf))

# Call
res <- glm(formula = ynorm ~ x, method = cirls.fit, Cmat = list(x = cpos),
  lb = list(x = 0), ub = list(x = Inf))
resdc <- do.call(glm, arglist)

# Check they are identical (save from call)
test_that("cirls can be used within do.call", {
  expect_equal(res[!names(res) %in% "call"], resdc[!names(res) %in% "call"])
})

# Check that the control argument has been modified in both cases
test_that("control argument is managed", {
  expect_null(res$control$Cmat)
  expect_null(resdc$control$Cmat)
})


#-----------------------
# Empty Cmat
#-----------------------

test_that("Everything works with empty Cmat", {

  # Check we have a warning when fitting with empty Cmat
  Cemp <- matrix(nrow = 0, ncol = p + 1)
  expect_warning(rese <-
      glm(formula = ynorm ~ x, method = cirls.fit, Cmat = Cemp))

  # Check we get same fit with unconstrained model
  resu <- glm(formula = ynorm ~ x)
  expect_equal(coef(rese), coef(resu))

  # Check inference runs fine
  expect_no_error(simulCoef(rese))
  expect_no_error(vcov(rese))
  expect_no_error(confint(rese))
})

