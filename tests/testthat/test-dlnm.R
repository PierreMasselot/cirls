################################################################################
#
# Tests for shape-constrained DLNM
#
################################################################################

library(dlnm)

#-----------------------------
# Simulate data
#-----------------------------

# Some parameters
maxlag <- 21

#----- Functions for simulation

# On variable dimension
fpara <- function(x) -.5 * sin(x * pi)

# Lag dimension
wdecay <- function(lag) exp(-lag/5)
wpeak <- function(lag) 5 * dnorm(lag, 5, 5)

# Full surface
fsurf <- function(x,lag) 0.1 * fpara(x) * wdecay(lag)

#----- Simulate data

# Use actual data for X
X <- london$tmean
X <- (X - min(X, na.rm = T)) / diff(range(X, na.rm = T))

# Compute cumulated effect
seqlag <- 0:maxlag
Q <- sapply(seqlag, function(l) c(rep(NA, l), X[1:(length(X) - l)]))
cumeff <- apply(Q, 1, function(hist) sum(fsurf(hist, seqlag)))

# Simulate response
set.seed(21)
suppressWarnings(y <- rpois(length(X), exp((log(20) + cumeff))))

#-----------------------------
# test shape-constrained DLNM
#-----------------------------

test_that("shape-constrained dlnm works", {

# List of type of basis functions
funlist <- c("strata", "bs", "ns")

# Loop over the types
for (vf in funlist) for (lf in funlist){
  # Define crossbasis
  cb <- crossbasis(X, lag = maxlag, argvar = list(fun = vf, df = 10),
    arglag = list(fun = lf, df = 5))

  # In what follows, dim is set to "var" by default (through `chkdim`)

  #----- Constraining the whole surface
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ shape(cb, shape = "inc"))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = .5)
  expect_true(all(diff(ccp$matfit) >= -sqrt(.Machine$double.eps)))

  #----- Only the overall
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ shape(cb, shape = "inc", overall = TRUE))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = .5)
  expect_true(all(diff(ccp$allfit) >= -sqrt(.Machine$double.eps)))

  #----- The overall but convex
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ shape(cb, shape = "cvx", overall = TRUE))

  # Prediction - each column (lag) should be non-decreasing
  at <- quantile(X, 0:10 / 10, na.rm = T)
  ccp <- crosspred(cb, cmodel, cen = .5, at = at)
  expect_true(all(diff(ccp$allfit, diff = 2) >= -sqrt(.Machine$double.eps)))

  #----- Slicing
  cmodel <- glm(y ~ cb, family = "quasipoisson", method = "cirls.fit",
    constr = ~ shape(cb, shape = "inc", slice = c(0, 5)))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = .5)
  expect_true(all(diff(ccp$matfit[,1]) >= -sqrt(.Machine$double.eps)))
  expect_true(all(diff(ccp$matfit[,5]) >= -sqrt(.Machine$double.eps)))
}

})


#-----------------------------
# test bound-constrained dlnm
#-----------------------------


test_that("bound constrained dlnm works", {

  # List of type of basis functions
funlist <- c("strata", "bs", "ns")

# Loop over the types
for (vf in funlist) for (lf in funlist){

  # Define crossbasis
  cb <- crossbasis(X, lag = maxlag, argvar = list(fun = vf, df = 10),
    arglag = list(fun = lf, df = 5))

  # In what follows, dim is set to "var" by default (through `chkdim`)

  #----- Constraining the whole surface
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ bound(cb))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = 0)
  expect_true(all(tail(ccp$matfit, 1) - 0 < sqrt(.Machine$double.eps)))

  #----- With another value
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ bound(cb, value = 1))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = 0)
  expect_true(all(tail(ccp$matfit, 1) - 1 < sqrt(.Machine$double.eps)))


  #----- Only the overall
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ bound(cb, overall = TRUE))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = 0)
  expect_equal(tail(ccp$allfit, 1), 0, ignore_attr = TRUE, tolerance = 10e-8)

  #----- Another value
  cmodel <- glm(y ~ cb, family = "quasipoisson",
    method = "cirls.fit", constr = ~ bound(cb, value = 1, overall = TRUE))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = 0)
  expect_equal(tail(ccp$allfit, 1), 1, ignore_attr = TRUE, tolerance = 10e-8)

  #----- Slicing
  cmodel <- glm(y ~ cb, family = "quasipoisson", method = "cirls.fit",
    constr = ~ bound(cb, slice = c(0, 5)))

  # Prediction - each column (lag) should be non-decreasing
  ccp <- crosspred(cb, cmodel, cen = 0)
  expect_equal(tail(ccp$matfit, 1)[1], 0, ignore_attr = TRUE, tolerance = 10e-8)
  expect_equal(tail(ccp$matfit, 1)[5], 0, ignore_attr = TRUE, tolerance = 10e-8)
}

})



#-----------------------------
# Constraing lags in DLNM
#-----------------------------

# # Define crossbasis
# cb <- crossbasis(X, lag = maxlag, argvar = list(fun = "bs", df = 10),
#   arglag = list(fun = "ns", df = 5))
#
# # Model
# cmodel <- glm(y ~ cb, family = "quasipoisson",
#     method = "cirls.fit", constr = ~ shape(cb, shape = "inc") +
#     shape(cb, shape = "dec", dim = "lag"))
# ccp <- crosspred(cb, cmodel, cen = 0)

# ################################################################################
# # DLMN test lags
#
# #--------------------------
# # Resimulate data
# #--------------------------
#
# # Create cumulative effect
# fsurf3 <- function(x,lag) 0.1 * fpara(x) * wpeak(lag)
# cumeff3 <- apply(Q, 1, function(hist) sum(fsurf3(hist, seqlag)))
#
# # Simulate response
# set.seed(31)
# suppressWarnings(y3 <- rpois(length(X), exp((log(20) + cumeff3))))
#
# #--------------------------
# # Models
# #--------------------------
#
# # Crossbasis
# cb <- crossbasis(X, lag = maxlag, argvar = list(fun = "bs", df = 10),
#     arglag = list(fun = "bs", df = 5))
#
# # Unconstrained model for comparison
# um <- glm(y3 ~ cb, family = "poisson")
# ucp <- crosspred(cb, um)
# plot(ucp, ptype = "slice", var = c(0, .5))
#
# #----- Constrained model
# cmdec <- glm(y3 ~ cb, family = "quasipoisson", method = "cirls.fit",
#     constr = ~ shape(cb, vshape = "dec", lshape = "dec"))
# cp <- crosspred(cb, cmdec)
# plot(cp, ptype = "slice", lag = c(0, 5, 10, 21), ci = "n")
# plot(cp, ptype = "overall")
# # coefmat <- matrix(coef(cmdec)[-1], ncol = 10)
#
#
# #----- Try to create a constraint matrix
# rng <- range(X, na.rm = T)
# dims <- attr(cb, "df")
# at <- seq(rng[1], rng[2], length.out = dims[1] + 2)
# at2 <- seq(rng[1], rng[2], length.out = 40)
# varbasis <- do.call(dlnm::onebasis, c(list(x = at2), attr(cb, "argvar")))
# cenbasis <- do.call(dlnm::onebasis, c(list(x = mean(rng)), attr(cb, "argvar")))
# varbasis <- scale(varbasis, center = cenbasis, scale = F)
# lagvec <- seq(0, 21, by = 1)
# lagbasis <- do.call(dlnm::onebasis, c(list(x = lagvec), attr(cb, "arglag")))
#
# Clag <- shapeConstr(lagbasis, shape = "dec")
# Cmat <- varbasis %x% Clag$Cmat
# checkCmat(Cmat, reduce = F)[1:2]
#
# cmtry <- glm(y3 ~ cb, family = "quasipoisson", method = "cirls.fit",
#     Cmat = list(cb = Cmat))
# Xpred <- dlnm:::mkXpred(type = "cb", basis = cb, at = at2, predvar = at2,
#   predlag = lagvec, cen = 0)
# matfit <- matrix(Xpred %*% coef(cmtry)[-1], length(at2), length(lagvec))
