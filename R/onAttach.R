################################################################################
# Startup message

.onAttach <- function(libname, pkgname){

  # get package info
  pkginfo <- utils::packageDescription("cirls")

  # Message
  msg <- paste0("This is cirls ", pkginfo$Version, ".\n",
    "For an introduction to the package see ?cirls.")
  packageStartupMessage(msg)
}
