# Simulate a simple dataset
set.seed(5)
x <- rnorm(100)
y <- x + rnorm(100)

# Run the model
mod <- glm(y ~ x, method = "cirls.fit", Cmat = list(x = 1), lb = 0)

# Loglikelihood and AIC
logLik(mod)
AIC(mod)

# AIC with a different degree of freedom
ll <- logLik(mod, df = "o")
AIC(ll)
