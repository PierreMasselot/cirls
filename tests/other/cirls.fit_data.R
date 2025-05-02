################################################################################
#
# Data to work on `cirls.fit`
#
################################################################################


# Number of obs
n <- 1000

# Coefficients
betas <- c(0, 1, 2, -1, 1)
p <- length(betas)

#----- Generate data

# Uniform values between 0 and 1
X <- matrix(rnorm(n * p), n, p)

# Linear predictor
eta <- 1 + X %*% betas

# Simulate responses
ynorm <- eta + rnorm(n, 0, .2)

# A constraint matrix
cinc <- diff(diag(p))

#----- Prepare everything for cirls.fit
mf <- model.frame(ynorm ~ X)
x <- model.matrix(mf)
y <- model.response(mf)
weights <- model.weights(mf)
start <- NULL
etastart <- NULL
mustart <- NULL
offset <- model.offset(mf)
family <- gaussian()
control <- list(Cmat = list(X = cinc))
intercept <- TRUE
singular.ok <- TRUE
mt <- terms(mf)
