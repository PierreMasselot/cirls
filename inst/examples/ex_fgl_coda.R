###############################################################################
# Compositional regression example with the fgl dataset

# Select the composition and use log
x <- as.matrix(fgl[,2:9])
z <- log(x)

# Fit model
model <- glm(RI ~ z, data = fgl, method = "cirls.fit",
  cons = ~ zerosum(z))

# Coefficients sum to 1
coef(model)
sum(coef(model)[-1])

# Confidence intervals
confint(model)

