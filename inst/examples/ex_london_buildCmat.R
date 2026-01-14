###############################################################################
# Examples of constraint specifications

#----- constr interface

# Fit an initial model and extract model.frame
mf <- model.frame(glm(death ~ pm10 + o3 + co, data = london))

# Simple non-negative constraint
buildCmat(mf, constr = ~ shape(pm10, "pos"))

# Several constraints
buildCmat(mf, constr = ~ shape(pm10, "pos") + shape(o3, "pos"))

# Different constraints
buildCmat(mf, constr = ~ shape(pm10, "pos") + zerosum(o3, co))

### In practice, this is done directly when calling glm
model <- glm(death ~ pm10 + o3 + co, data = london, family = "quasipoisson",
  method = "cirls.fit", constr = ~ shape(pm10, "pos"))
model[c("Cmat", "lb", "ub")]
buildCmat(mf, constr = ~ shape(pm10, "pos"))

#----- Specific terms

# Simple bound constraint
buildCmat(mf, lb = list("pm10" = 0))

# Equivalent of zerosum
buildCmat(mf, Cmat = list("o3;co" = matrix(1, 1, 2)), ub = list("o3;co" = 0))

### In practice, this is done directly when calling glm
model <- glm(death ~ pm10 + o3 + co, data = london, family = "quasipoisson",
  method = "cirls.fit", lb = list("pm10" = 0))
model[c("Cmat", "lb", "ub")]
buildCmat(mf, lb = list("pm10" = 0))

#----- Both options simultaneously

# Bound constraint and zerosum
buildCmat(mf, lb = list("pm10" = 0), constr = ~  zerosum(o3, co))

# Same as before, should be done in glm
model <- glm(death ~ pm10 + o3 + co, data = london, family = "quasipoisson",
  method = "cirls.fit",
  lb = list("pm10" = 0), constr = ~  zerosum(o3, co))
model[c("Cmat", "lb", "ub")]

#----- Fully specified matrix

# Simple bound constraints
Cmat = cbind(0, diag(2), 0)
buildCmat(mf, Cmat = Cmat)

# Zerosum constraint
Cmat <- t(c(0, rep(1, 3)))
buildCmat(mf, Cmat = Cmat, lb = 0, ub = 0)

#----- Unconstrained model
buildCmat(mf)
