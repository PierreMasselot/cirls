#------------------------------
# Generate some data
#------------------------------

#----- Parameters

# Number of obs
n <- 1000

# Coefficients
betas <- c(1, 1)
p <- length(betas)

#----- Generate data

# Uniform values between 0 and 1
x <- matrix(rnorm(n * p), n, p)

# Linear predictor
eta <- x %*% betas

# Simulate responses
y <- eta + rnorm(n, 0, .2)

#----- Fit model
res <- glm(y ~ 0 + x, method = cirls.fit, Cmat = list(x = diff(diag(p))))


#------------------------------
# Simulate
#------------------------------

#----- With constraints

# Simulate
simcons <- simulCoef(res, nsim = 1000, seed = 5)

# Plot
plot(simcons)
abline(a = 0, b = 1, col = 2)

#----- Unconstrained

# Simulate
simun <- simulCoef(res, nsim = 1000, seed = 5, cons = F)

# Plot
plot(simun)
abline(a = 0, b = 1, col = 2)
