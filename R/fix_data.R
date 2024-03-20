## changes to blavaan-exported Stan data, for SAM
fix_data <- function(blavobject, standata) {
  freemats <- lavInspect(blavobject, "free")
  lavpartable <- parTable(blavobject)

  ## pick up blocks of psi, including "blocks" of one entry:
  if (inherits(freemats[[1]], "matrix")) freemats <- list(freemats)

  if ("psi" %in% names(freemats[[1]])) {
    blkinfo <- lapply(freemats, function(x) blkdiag(x$psi, attributes(freemats)$header))
    blkse <- do.call("rbind", lapply(blkinfo, function(x) x$blkse))

    if (nrow(blkse) > 0) {
      blkse <- blkse[(blkse[,3] == 1), , drop = FALSE] # for joint model: & (blkse[,2] - blkse[,1] > 1), , drop = FALSE]
      blksizes <- blkse[,2] - blkse[,1] + 1
      ublksizes <- unique(blksizes)
      ublksizes <- ublksizes[order(ublksizes)]

      ## for now, we only support 5 unique block dimensions because each dimension requires
      ## a separate parameter specification in Stan
      if (length(ublksizes) > 5 | length(ublksizes) == 0) {
        blkinfo <- NULL
        blkpsi <- FALSE
        standata$nblk <- array(0, dim = 5)
        standata$psidims <- array(3, dim = 5)
        standata$blkse <- matrix(nrow = 0, ncol = 7)
      } else {
        blkgrp <- rep(1:length(blkinfo), times = sapply(blkinfo, function(x) nrow(x$blkse)))
        arrayidx <- as.numeric(as.factor(ublksizes))
        dupsiz <- duplicated(blksizes)
        blkidx <- rep(NA, nrow(blkse))
        for (i in 1:length(ublksizes)) {
          sizeidx <- blksizes == ublksizes[i]
          blkidx[sizeidx] <- cumsum(dupsiz[sizeidx]) + 1
        }
        blkse <- cbind(blkse, blkgrp, arrayidx, blkidx, rep(1, nrow(blkse))) ## col 7 is for priors, handled later

        nblk <- c(summary(factor(blksizes)), rep(0, 5 - length(ublksizes)))
        standata$nblk <- array(nblk, dim = 5)
        psidims <- c(ublksizes, rep(3, 5 - length(ublksizes)))
        standata$psidims <- array(psidims, dim = 5)
        standata$blkse <- blkse
      }
      
      ## for correlations in blocks, replace beta with lkj and deal with lkj priors
      for (b in 1:nrow(blkse)) {
        lkjrows <- with(lavpartable, which(group == blkse[b, 'blkgrp'] & mat == "lvrho" &
                                           row >= blkse[b,1] & col <= blkse[b,2] & row != col))
        nolkj <- lkjrows[!grepl("lkj", lavpartable$prior[lkjrows])]
        if (length(nolkj) > 0) {
          lavpartable$prior[nolkj] <- gsub("(\\w+)\\(([^,]+),([^)]+)\\)", "lkj_corr(\\2)", lavpartable$prior[nolkj])
        }
        lptrow <- with(lavpartable, which(row == blkse[b,1] & col == (blkse[b,1] + 1) &
                                          group == blkse[b,4] & mat == "lvrho"))
        if (length(lptrow) > 0) {
          blkse[b,7] <- as.numeric(gsub("(\\w+)\\(([^,]+)\\)", "\\2", lavpartable$prior[lptrow]))
        } else {
          blkse[b,7] <- 1
        }
        lavpartable$prior[lkjrows] <- lavpartable$prior[lptrow]
      }
    }
  }

  ## add extra initial values
  inits <- blavobject@external$inits
  for (i in 1:length(inits)) {
    inits[[i]]$Psi_sd_tmp <- inits[[i]]$Psi_sd_free
  }

  list(standata = standata, lpt = lavpartable, inits = inits)
}
    
