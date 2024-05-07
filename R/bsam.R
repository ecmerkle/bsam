## BSAM: Bayesian Structural After Measurement
## uses a new Stan file, then uses blavaan for everything else

bsam <- function(...) {

  dotdotdot <- list(...)

  if (!("mcmcextra" %in% names(dotdotdot))) {
    mcmcextra <- list(dosam = TRUE)
  } else {
    mcmcextra <- dotdotdot$mcmcextra
    mcmcextra$dosam <- TRUE
  }

  mc <- match.call()
  mc[[1L]] <- quote(bsem)
  mc$mcmcextra <- mcmcextra

  eval(mc, parent.frame())

}
