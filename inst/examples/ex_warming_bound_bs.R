###############################################################################
# Bound constraints for bs: warming dataset

library(splines)

# Force the fit to reach the last value
# Remove the intercept, so that it is included in the constraint
basis <- bs(warming$year, df = 10)
mod <- glm(anomaly ~ 0 + basis, data = warming, method = "cirls.fit",
  constr = ~ bound(basis, value = 0.75))

# Plot result
plot(anomaly ~ year, data = warming, xlab = "", ylab = "Temperature anomaly")
lines(warming$year, predict(mod), col = 3, lwd = 2)
