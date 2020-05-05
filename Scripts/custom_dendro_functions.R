brf_RWB <- function (g, p, sb, vnames, ci = 0.05, NN = NNN) 
{
  n <- length(g)
  m <- dim(p)[2]
  param.matrix <- matrix(NA, nrow = m, ncol = NN)
  if (sb) {
    pb <- txtProgressBar(min = 1, max = NN, style = 3)
  }
  for (i in 1:NN) {
    boot.sample <- sample(1:n, n, replace = TRUE)
    boot.g <- g[boot.sample]
    boot.p <- p[boot.sample, ]
    boot.g <- (boot.g - mean(boot.g))/sd(boot.g)
    boot.p <- apply(boot.p, 2, function(x) {
      (x - mean(x))/sd(x)
    })
    cor.mat <- cor(boot.p)
    eigen.decomp <- eigen(cor.mat)
    eigenvectors <- eigen.decomp$vectors
    eigenvalues <- eigen.decomp$values
    cumprods <- cumprod(eigenvalues)
    reduced.eigenvectors <- eigenvectors[, cumprods > 1]
    pc.scores <- boot.p %*% reduced.eigenvectors
    k <- qr.solve(pc.scores, boot.g)
    zeros <- rep(0, length(which(cumprods < 1)))
    k <- c(k, zeros)
    b <- eigenvectors %*% k
    param.matrix[, i] <- b
    if (sb) 
      setTxtProgressBar(pb, i)
  }
  brf.coef <- apply(param.matrix, 1, median)
  if (ci == 0.05) {
    ci.lower <- apply(param.matrix, 1, function(x) {
      sort(x)[NN*0.05]
    })
    ci.upper <- apply(param.matrix, 1, function(x) {
      sort(x)[NN*0.95]
    })
  }
  else {
    if (ci == 0.01) {
      ci.lower <- apply(param.matrix, 1, function(x) {
        sort(x)[NN*0.01]
      })
      ci.upper <- apply(param.matrix, 1, function(x) {
        sort(x)[NN*0.99]
      })
    }
    else {
      if (ci == 0.1) {
        ci.lower <- apply(param.matrix, 1, function(x) {
          sort(x)[NN*0.1]
        })
        ci.upper <- apply(param.matrix, 1, function(x) {
          sort(x)[NN*0.9]
        })
      }
      else {
        stop("`ci` must be either 0.1, 0.05, or 0.01.")
      }
    }
  }
  is.sig <- logical(m)
  for (i in 1:m) {
    if (sign(ci.upper[i]) != sign(ci.lower[i])) {
      is.sig[i] <- FALSE
    }
    else {
      if (abs(brf.coef[i]) > abs((abs(ci.upper[i]) - abs(ci.lower[i]))/2)) {
        is.sig[i] <- TRUE
      }
      else {
        is.sig[i] <- FALSE
      }
    }
  }
  out <- data.frame(coef = brf.coef, significant = is.sig, 
                    ci.lower = ci.lower, ci.upper = ci.upper)
  rownames(out) <- colnames(p)
  if (sb) 
    close(pb)
  attributes(out)$npar <- attributes(p)$npar
  attributes(out)$vnames <- vnames
  out
}

dcc_RWB <- function (chrono, clim, method = "response", start = -6, end = 9, 
                     timespan = NULL, vnames = NULL, sb = TRUE, boot = TRUE, ci = 0.05, NNN = 10000) 
{
  month.ids <- c(-1:-12, 1:12)
  errormsg1 <- "start and end have to define an interval in [-1, -2, ..., -12, 1, 2, ..., 12]."
  if (!is.element(start, month.ids) || !is.element(end, month.ids) || 
      which(month.ids == start) > which(month.ids == end)) {
    stop(errormsg1)
  }
  clim <- climdispatch(clim)
  chrono.years <- as.numeric(row.names(chrono))
  clim.years <- sort(unique(clim[, 1]))
  if (chrono.years[1] <= clim.years[1]) {
    overlap <- na.omit(clim.years[match(chrono.years, clim.years)])
  }
  else {
    overlap <- na.omit(chrono.years[match(clim.years, chrono.years)])
  }
  if (is.null(timespan)) {
    start.year <- overlap[1]
    end.year <- tail(overlap, 1)
  }
  else {
    if (start > 0) {
      if (!is.element(timespan[1], overlap) || !is.element(timespan[2], 
                                                           overlap)) {
        errormsg3 <- paste("timespan has to be between ", 
                           overlap[1], " and ", tail(overlap, 1), " for start dates in current year.", 
                           sep = "")
        stop(errormsg3)
      }
      else {
        start.year <- timespan[1]
        end.year <- timespan[2]
      }
    }
    else {
      if (!is.element(timespan[1], overlap) || !is.element(timespan[2], 
                                                           overlap)) {
        errormsg4 <- paste("timespan has to be between ", 
                           overlap[1] + 1, " and ", tail(overlap, 1), 
                           " for start dates in previous year.", sep = "")
        stop(errormsg4)
      }
      else {
        start.year <- timespan[1]
        end.year <- timespan[2]
      }
    }
  }
  if (start < 0 && is.na(match((start.year - 1), clim.years))) {
    offset <- 1
  }
  else {
    offset <- 0
  }
  if (start < 0) {
    interval.clim <- (start.year - 1 + offset):end.year
    interval.chrono <- (start.year + offset):end.year
  }
  else {
    interval.clim <- (start.year + offset):end.year
    interval.chrono <- (start.year + 1 + offset):end.year
  }
  if (start * end > 0) {
    no.params <- (dim(clim)[2] - 2) * length(start:end)
  }
  else {
    no.params <- (dim(clim)[2] - 2) * length(start:end) - 
      1
  }
  overlap.size <- length(start.year:end.year)
  if (no.params > overlap.size) {
    win.size.msg <- paste("Overlapping time span of chrono and climate records is smaller than number of parameters! Consider adapting the number of parameters to a maximum of ", 
                          overlap.size, ".", sep = "")
    stop(win.size.msg)
  }
  a <- as.numeric(rownames(chrono)) %in% interval.chrono
  b <- clim[, 1] %in% interval.clim
  chrono.trunc <- chrono[a, 1]
  clim.trunc <- clim[b, ]
  p <- pmat(clim.trunc, start, end, vnames)
  METHOD <- match.arg(method, c("response", "correlation"))
  if (METHOD == "response") {
    if (boot) {
      dc <- brf_RWB(chrono.trunc, p, sb = sb, vnames = vnames, NN = NNN)
    }
    else {
      dc <- nbrf(chrono.trunc, p, vnames = vnames)
    }
  }
  if (METHOD == "correlation") {
    if (boot) {
      dc <- bcf(chrono.trunc, p, sb = sb, vnames = vnames, 
                ci = ci)
    }
    else {
      dc <- nbcf(chrono.trunc, p, vnames = vnames)
    }
  }
  cat("time span considered:", start.year, "-", end.year, "\n")
  dc
}

