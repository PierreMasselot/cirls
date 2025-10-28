################################################################################
#
#  Test what happens when model multicollinear
#
################################################################################

# Generate perfectly colinear variables
x1 <- rnorm( 100 )
x2 <- 2 * x1
y <- rnorm( 100 )

# Check that model fitting returns warning
test_that("multicollinear model works", {
  expect_warning(
    mod <- glm(y ~ x1 + x2, method = cirls.fit, Cmat = cbind(0, diag(2)),
      lb = 0, ub = Inf),
    "constraints removed")
  expect_lt(length(na.omit(coef(mod))), 3)

  # Vcov
  v <- vcov(mod)
  aliased <- summary(mod)$aliased
  expect_true(all(is.na(
    v[row(v) %in% which(aliased) | col(v) %in% which(aliased)])))
  expect_false(any(is.na(
    v[!(row(v) %in% which(aliased) | col(v) %in% which(aliased))])))

  # Confint
  ci <- confint(mod)
  expect_lt(nrow(na.omit(ci)), 3)
  expect_true(all(is.na(ci[aliased,])))
})
