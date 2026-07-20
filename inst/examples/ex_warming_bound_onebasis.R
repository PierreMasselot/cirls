###############################################################################
# Bound constraints for onebasis: warming dataset

library(dlnm)

# Force the fit to reach the last value
# Remove the intercept, so that it is included in the constraint
basis <- onebasis(warming$year, fun = "strata", df = 10, intercept = TRUE)
mod <- glm(anomaly ~ 0 + basis, data = warming, method = "cirls.fit",
  constr = ~ bound(basis, value = 0.75))

# Plot result
plot(anomaly ~ year, data = warming, xlab = "", ylab = "Temperature anomaly")
lines(warming$year, predict(mod), col = 2, lwd = 2)
