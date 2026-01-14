###############################################################################
# Compositional regression example with the fgl dataset

# Add a small increment to zeros
x <- as.matrix(fgl[,2:9])
x[x == 0] <- min(x[x > 0]) / 2

# Normalise to unit sum and take log
z <- log(x / rowSums(x))

# Fit model
model <- glm(RI ~ z, data = fgl, method = "cirls.fit",
  cons = ~ zerosum(z))

# Coefficients and confidence intervals
coef(model)
confint(model)
