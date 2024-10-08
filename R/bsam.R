## BSAM: Bayesian Structural After Measurement
## uses a new Stan file, then uses blavaan for everything else

bsam <- function(..., ngibbs = 50L) {

  mc <- match.call()
  mc[[1L]] <- quote(bsem)
  
  dotdotdot <- list(...)
  if (!("mcmcextra" %in% names(dotdotdot))) {
    mcmcextra <- list(dosam = TRUE, data = list(ngibbs = ngibbs, fullpsi_c = 0), monitor = "PS")
  } else {
    mcmcextra <- dotdotdot$mcmcextra
    mcmcextra$dosam <- TRUE
    mcmcextra$monitor <- c(mcmcextra$monitor, "PS")
    mcmcextra$data <- c(mcmcextra$data, list(ngibbs = ngibbs, fullpsi_c = 0L))
  }
  mc$mcmcextra <- mcmcextra
  
  mc$meanstructure <- TRUE
  if ("meanstructure" %in% names(dotdotdot)) {
    if (!dotdotdot$meanstructure) {
      warning("bsam WARNING: meanstructure has no effect and will be set to TRUE")
    }
  }

  eval(mc, parent.frame())
}
