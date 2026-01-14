###############################################################################
# Example of AIC and logLik

#----- Fit warming model

# Define decadal strata (as a factor)
warming$decade <- factor(10 * floor(warming$year / 10))

# Non-decreasing constraint on decadal strata
model <- glm(anomaly ~ decade, data = warming, method = "cirls.fit",
  cons = ~ shape(decade, "inc"))

#----- Log-Likelihood and information criteria

# Extract log-likelihood with different degrees of freedom
lle <- logLik(model)
llo <- logLik(model, df = "o")
llu <- logLik(model, df = "u")

# The AIC function can be used directly on the model: uses default df
# NB: the small differences come from the simulations for edf computation
AIC(model)
AIC(lle)

# To compute AIC with different df, need to go through logLik
AIC(llo)
AIC(llu)
