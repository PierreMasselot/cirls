################################################################################
#
# Creates the first derivative matrix of B-splines with custom knots and order
#
################################################################################

dm <- function(l, knots, ord){

  # Compute lags between knots and invert
  kind <- (l + 1):(length(knots) - l)
  dlknots <- diff(knots[kind], ord - l)

  # Average knot length (to avois very small or large numbers in Cmat)
  avkl <- mean(diff(knots))

  # Compute the difference matrix
  diff(diag(length(dlknots) + 1)) / dlknots * (ord - l) * avkl
}
