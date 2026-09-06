###############################################################################
# Monotone strata levels with the warming dataset

### Typical usage

# Non-decreasing constraint on decadal strata
model <- glm(anomaly ~ decade, data = warming, method = "cirls.fit",
  cons = ~ shape(decade, "inc"))

# Extract covariance matrix and confidence intervals
v <- vcov(model)
ci <- confint(model)

### Using the sim.cirls objects

# Simulate coefficients directly
sims <- simulCoef(model, nsim = 1000, seed = 4)

# Get covariance matrix and confidence intervals
vs <- vcov(sims)
cis <- confint(sims)

# Check that we get the same results as the typical usage
v2 <- vcov(model, seed = 4)
ci2 <- confint(model, seed = 4)
identical(vs, v2)
identical(cis, ci2)
