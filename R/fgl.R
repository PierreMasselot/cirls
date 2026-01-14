################################################################################
# Dataset documentation: fgl

#' Measurement of forensic glass fragments
#'
#' Chemical composition of fragments of glass collected in forensic work.
#'
#' @format A data frame of 214 observations and 10 variables. Apart from `RI` and `type` all other variables represent a percentage of the given chemical component.
#' \describe{
#'  \item{RI}{Refractive index; more precisely the refractive index is 1.518xxxx.}
#'  \item{Na}{Sodium}
#'  \item{Mg}{Manganese}
#'  \item{Al}{Aluminium}
#'  \item{Si}{Silicon}
#'  \item{K}{Potassium}
#'  \item{Ca}{Calcium}
#'  \item{Ba}{Barium}
#'  \item{Fe}{Iron}
#'  \item{type}{The fragments were originally classed into seven types, one of which was absent in this dataset. The categories which occur are window float glass (`WinF`: 70), window non-float glass (`WinNF`: 76), vehicle window glass (`Veh`: 17), containers (`Con`: 13), tableware (`Tabl`: 9) and vehicle headlamps (`Head`: 29)}
#' }
#'
#' @source
#' Dataset exctracted from the `MASS` package.
#'
#' Venables, W. N. and Ripley, B. D. (2002) *Modern Applied Statistics with S*. Fourth edition. Springer.
#'
#' @example inst/examples/ex_fgl_coda.R
#'
"fgl"

