################################################################################
# Dataset documentation: london

#' Daily mortality, temperature and pollution data in London
#'
#' This dataset contains daily time series of all-cause mortality for several age groups, temperature (maximum, mean and minimum) and air pollutants (PM10, ozone and carbon monoxide) for the period of 1996-2003.
#'
#' @format A data frame with 2922 observations and 12 variables:
#' \describe{
#'  \item{date}{`Date` variable}
#'  \item{death}{Counts of all cause and all-age mortality}
#'  \item{age0_64}{Counts of mortality for age group 0-64}
#'  \item{age65_74}{Counts of mortality for age group 65-74}
#'  \item{age75_84}{Counts of mortality for age group 75-84}
#'  \item{age85plus}{Counts of mortality for age group 85+}
#'  \item{tmax}{Daily maximum temperature (in Celsius degree)}
#'  \item{tmean}{Daily mean temperature (in Celsius degree)}
#'  \item{tmin}{Daily minimum temperature (in Celsius degree)}
#'  \item{pm10}{Daily mean coarse particulate matter concentration (\eqn{\mu g/m^3})}
#'  \item{o3}{Daily mean ozone concentration (\eqn{\mu g/m^3})}
#'  \item{co}{Daily mean carbon monoxide concentration (\eqn{\mu g/m^3})}
#' }
#'
#' @source
#' Office for National Statistics
#'
#' British Atmospheric Data Centre
#'
#' @example inst/examples/ex_london_nonneg.R
#'
"london"
