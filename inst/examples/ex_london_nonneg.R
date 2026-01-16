###############################################################################
# Non-negative coefficient example with the london dataset

library(splines)

### Association between CO and mortality
# Model includes spline of date and day-of-week
# Nonnegative constraint on CO
model <- glm(death ~ co + weekdays(date) + ns(date, df = 7*8),
  data = london, family = "quasipoisson",
  method = "cirls.fit", cons = ~ shape(co, "pos"))

# Coefficient and confidence interval
coef(model)[2]
confint(model)[2,]

# Comparing to an unconstrained model
umodel <- uncons(model)
coef(umodel)[2]
confint(umodel)[2,]
