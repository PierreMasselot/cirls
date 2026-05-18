###############################################################################
# Constrained DLNM: var dimension

\dontrun{

library(dlnm)
library(splines)

#----------------
# Association between temperature and mortality
#----------------

# Number of years and dow
ny <- length(unique(format(london$date, "%Y")))

# Create a flexible crossbasis
cb <- crossbasis(london$tmean, lag = 21,
  argvar = list(fun = "bs", degree = 2, df = 10),
  arglag = list(fun = "ns", knots = logknots(21, df = 5)))

#----- Fit several models

# In the following, `dim` is not shown as it is set to "var" by default

# A nondecreasing model
icm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "inc"))

# A convex model
ccm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "cvx"))

# A convex model with a constraint on overall only
ocm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "cvx", overall = TRUE))

# A convex model effectively constraining only lags 0 and 1
subcm <- glm(age0_64 ~ cb + ns(date, df = 7 * ny) + dow, data = london,
  family = "quasipoisson", method = "cirls.fit",
  constr = ~ shape(cb, shape = "cvx", slice = c(0, 1)))

#----- Plot results

# Plot non-decreasing
icp <- crosspred(cb, icm, cen = 20)
par(mfrow = c(1, 2))
plot(icp, main = "Non-decreasing")
plot(icp, ptype = "overall")

# Plot convex
ccp <- crosspred(cb, ccm, cen = 20)
par(mfrow = c(1, 2))
plot(ccp, main = "Convex")
plot(ccp, ptype = "overall")

# This one is restricted to overall but not for every single lag
ocp <- crosspred(cb, ocm, cen = 20)
par(mfrow = c(1, 2))
plot(ocp, main = "Overall convex")
plot(ocp, ptype = "overall")

# Restricted to lags 0 and 1 only
subcp <- crosspred(cb, subcm, cen = 20)
par(mfrow = c(1, 2))
plot(subcp, main = "Convex")
plot(subcp, ptype = "overall")

}
