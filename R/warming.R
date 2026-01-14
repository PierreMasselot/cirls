################################################################################
# Dataset documentation: warming

#' Global temperature anomaly
#'
#' Provides annual global temperature anomaly compared to the reference period 1961-1990 from 1850 to 2015
#'
#' @format
#' A data frame of 166 observations (years) and two? variables:
#' \describe{
#'  \item{year}{The year.}
#'  \item{anomaly}{The difference between the global annual temperature and the average over 1961-1990.}
#' }
#'
#' @source
#' Jones, P.D., Parker, D.E., Osborn, T.J., Briffa, K.R., 2000. Global and Hemispheric Temperature Anomalies: Land and Marine Instrumental Records (1850 - 2015). [DOI:10.3334/CDIAC/CLI.002](https://doi.org/10.3334/CDIAC/CLI.002)
#'
#' @example inst/examples/ex_warming_factor.R
#'
"warming"
