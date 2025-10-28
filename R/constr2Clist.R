################################################################################
#
# Create list of constraint matrices from constr formula
#
################################################################################

constr2Clist <- function(var, label, int, mf) {
  
  # prepare call
  ctcstr <- paste0(var[[1]], "Constr")
  if(!exists(ctcstr)) 
    stop(sprintf("constraint method for %s does not exist", label))
  var[[1]] <- str2lang(ctcstr)
  if(int) var$intercept <- TRUE
  
  # Evaluate and return
  eval(var, mf)
}
