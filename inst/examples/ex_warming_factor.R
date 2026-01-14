###############################################################################
# Monotone strata levels with the warming dataset

### Fit the model

# Define decadal strata (as a factor)
warming$decade <- factor(10 * floor(warming$year / 10))

# Non-decreasing constraint on decadal strata
model <- glm(anomaly ~ decade, data = warming, method = "cirls.fit",
  cons = ~ shape(decade, "inc"))

# Plot result
plot(anomaly ~ year, data = warming, xlab = "", ylab = "Temperature anomaly")
lines(warming$year, predict(model), col = 2, lwd = 2)

### Extract results

# Coefficients and confidence intervals
betas <- coef(model)
v <- vcov(model)
cis <- confint(model)

# Degrees of freedom: represent number of strata change
# ?edf
edf(model)

### Compare with an unconstrained model

# Refit the model
umodel <- uncons(model)

# Add result
lines(warming$year, predict(umodel), col = 3, lwd = 2, lty = 2)

