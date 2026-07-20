###############################################################################
# Bound constraints example: warming dataset

library(splines)

# Force the fit to reach the last value
# Remove the intercept, so that it is included in the constraint
model <- glm(anomaly ~ 0 + decade, data = warming, method = "cirls.fit",
  constr = ~ bound(decade, value = 0.75))

# Same but with splines
basis <- ns(warming$year, df = 10)
splmodel <- glm(anomaly ~ 0 + basis, data = warming, method = "cirls.fit",
  constr = ~ bound(basis, value = 0.75))

# Plot result
plot(anomaly ~ year, data = warming, xlab = "", ylab = "Temperature anomaly")
lines(warming$year, predict(model), col = 2, lwd = 2)
lines(warming$year, predict(splmodel), col = 3, lwd = 2)
