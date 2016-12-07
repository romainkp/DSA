.onLoad <- function(...) {
  ## this is for the model fitting routines; really internal debugging
  set.DSA.fit.debug(10) 

  ## this is for the DSA itself, often useful to trace.
  setDSAMessageLevel(-1)

  ## this is for the comparison.
  .C("DSA_PACK_setEpsilonCompare", as.double(.Machine$double.eps*100))
}

.onAttach <- function(...) {
  library(modelUtils)
  cat("Welcome to the DSA \n")

  if(hasDSASubversionInfo()) {
    getDSASubversionInfo()
  }
  else {
    print.package.message(NULL)
  }
}

.onUnload <- function(...) {
  library.dynam.unload("DSA", libpath)
}
