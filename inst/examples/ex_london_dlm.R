###############################################################################
# Constrained distributed lag models (DLM)

library(dlnm)
library(splines)

#----------------
# Association between PM10 and mortality
#----------------

# Number of years and dow
ny <- length(unique(format(london$date, "%Y")))

# Create a flexible crossbasis
dlm <- crossbasis(london$pm10, lag = 15, argvar = list(fun = "lin"),
  arglag = list(fun = "ns", knots = logknots(15, df = 4)))

#----- Fit models

# We always have to set `dim = "lag"` because it is set to "var" by default

# Fit with bound constraint, value = 0 by default
bm <- glm(age0_64 ~ dlm + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ bound(dlm, dim = "lag"))

# Get a decaying lag-response function with a decreasing constraint
dm <- glm(age0_64 ~ dlm + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(dlm, dim = "lag", shape = "dec"))

# Both can be used at the same time
dbm <- glm(age0_64 ~ dlm + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(dlm, dim = "lag", shape = "dec") + bound(dlm, dim = "lag"))

#----- Plot to compare

# Compare to unconstrained model
um <- glm(age0_64 ~ dlm + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson")

# Use crosspred to get the lag-response function
ucp <- crosspred(dlm, um, cen = 0, at = 10)
bcp <- crosspred(dlm, bm, cen = 0, at = 10)
dcp <- crosspred(dlm, dm, cen = 0, at = 10)
dbcp <- crosspred(dlm, dbm, cen = 0, at = 10)

# Now plot
par(mfrow = c(2,2))
plot(ucp, ptype = "slices", lwd = 2, var = 10, main = "Unconstrained")
plot(bcp, ptype = "slices", lwd = 2, var = 10, col = 2,
  main = "Bound")
plot(dcp, ptype = "slices", lwd = 2, var = 10, col = 3,
  main = "Decreasing")
plot(dbcp, ptype = "slices", lwd = 2, var = 10, col = 4,
  main = "Bound & decreasing")


