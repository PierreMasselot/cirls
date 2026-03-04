###############################################################################
# Full isotonic regression with warming dataset

### Fit the model

# Factor for year
warming$fyear <- factor(warming$year)

# Non-decreasing constraint on decadal strata
model <- glm(anomaly ~ fyear, data = warming, method = "cirls.fit",
  cons = ~ shape(fyear, "inc"))
