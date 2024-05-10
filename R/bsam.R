## BSAM: Bayesian Structural After Measurement
## uses a new Stan file, then uses blavaan for everything else

bsam <- function(..., ngibbs = 50L) {

  dotdotdot <- list(...)

  if (!("mcmcextra" %in% names(dotdotdot))) {
    mcmcextra <- list(dosam = TRUE, data = list(ngibbs = ngibbs, fullpsi_c = 0), monitor = "PS")
  } else {
    mcmcextra <- dotdotdot$mcmcextra
    mcmcextra$dosam <- TRUE
    mcmcextra$monitor <- c(mcmcextra$monitor, "PS")
    mcmcextra$data <- c(mcmcextra$data, list(ngibbs = ngibbs, fullpsi_c = 0L))
  }

  mc <- match.call()
  mc[[1L]] <- quote(bsem)
  mc$mcmcextra <- mcmcextra

  eval(mc, parent.frame())

}
