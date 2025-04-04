# Simulate a simple dataset
set.seed(5)
x <- rnorm(100)
y <- x + rnorm(100)

#### Model with a sensible constraint
# Reduces edf compared to udf as the constraint is sometimes active
mod1 <- glm(y ~ x, method = "cirls.fit", Cmat = list(x = 1), lb = 1)
edf(mod1)

#### Model with an almost surely binding constraint
# In this case edf is very close to odf as the constraint is often active
mod2 <- glm(y ~ x, method = "cirls.fit", Cmat = list(x = 1), lb = 1.5)
edf(mod2)

#### Model with an irrelevant constraint
# Here the constraint is useless and edf is equal to unconstrained df
mod3 <- glm(y ~ x, method = "cirls.fit", Cmat = list(x = 1), lb = -5)
edf(mod3)
