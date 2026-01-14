##############################################################################
# Monotone splines with the warming dataset

library(splines)

# Define a natural spline basis
basis <- ns(warming$year, df = 10)

### Constrained model

# Non-decreasing splines
model <- glm(anomaly ~ basis, data = warming, method = "cirls.fit",
  cons = ~ shape(basis, "inc"))

# Plot result
plot(anomaly ~ year, data = warming, xlab = "", ylab = "Temperature anomaly")
lines(warming$year, predict(model), col = 2, lwd = 2)

### Unconstrained model
umodel <- uncons(model)
lines(warming$year, predict(umodel), col = 3, lwd = 2, lty = 2)

