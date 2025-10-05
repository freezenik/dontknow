library("mvtnorm")

## Helper functions.
invrhogit <- function(eta) {
  eta / sqrt(1 + eta^2)
}

rhogit <- function(mu) {
  mu / sqrt(1 - mu^2)
}

cut2inc <- function(cut) {
  inc <- cut
  if(length(cut) > 1L) {
    inc[2:length(cut)] <- log(diff(cut))
  }
  return(inc)  
}

inc2cut <- function(inc) {
  cut <- inc
  if(length(inc) > 1L) {
    cut[2L:length(inc)] <- cut[1] + cumsum(exp(inc[2L:length(inc)]))
  }
  return(cut)  
}

## Log-likelihood function.
loglik <- function(eta1, eta2, rho, alpha1, alpha2, y, log = TRUE) {
  y1 <- y[, 1L]
  y2 <- y[, 2L]

  n  <- length(y1)
  ll <- numeric(n)

  rho <- pmin(0.999999, pmax(-0.999999, rho[1L]))

  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  chSigma <- chol(Sigma)
  chSigma <- mvtnorm::as.ltMatrices(t(chSigma))

  idx1 <- which(y1 == 1L)
  if(length(idx1)) {
    ll[idx1] <- pnorm(alpha1[2] - eta1[idx1], lower.tail = FALSE, log.p = TRUE)
  }

  idx0 <- which(y1 == 0L)
  if(length(idx0)) {
    lower <- rbind(rep(-Inf, length(idx0)),
                   alpha2[y2[idx0] + 1L] - eta2[idx0])
    upper <- rbind(alpha1[2] - eta1[idx0],
                   alpha2[y2[idx0] + 2L] - eta2[idx0])

    ll[idx0] <- lpmvnorm(lower = lower, upper = upper,
                         mean = c(0, 0), chol = chSigma,
                         logLik = FALSE, M = 250)
  }

  if(!log)
    ll <- exp(ll)

  ll
}


## Families.
dk3 <- function(...) {
  f <- list(
    "family" = "dk3",
    "names" = c("mu1", "mu2", "rho", "cutinc11", "cutinc21", "cutinc22"),
    "links" = c("identity", "identity", "rhogit", "identity", "identity", "identity"),
    "d" = function(y, par, log = FALSE) {
      alpha1 <- c(-Inf, inc2cut(par$cutinc11[1]), Inf)
      alpha2 <- c(-Inf, inc2cut(c(par$cutinc21[1], par$cutinc22[1])), Inf)
      
      eta1 <- par$mu1
      eta2 <- par$mu2
      
      rho <- par$rho
      
      d <- loglik(eta1 = eta1, eta2 = eta2, rho = rho, alpha1 = alpha1,
        alpha2 = alpha2, y = y, log = log)
      
      return(d)
    }
  )
  class(f) <- "gamlss2.family"
  return(f)
}

dk4 <- function(...) {
  f <- list(
    "family" = "dk4",
    "names" = c("mu1", "mu2", "rho", "cutinc11", "cutinc21", "cutinc22", "cutinc23"),
    "links" = c("identity", "identity", "rhogit", "identity",
      "identity", "identity", "identity"),
    "d" = function(y, par, log = FALSE) {
      alpha1 <- c(-Inf, inc2cut(par$cutinc11[1]), Inf)
      alpha2 <- c(-Inf,
        inc2cut(c(par$cutinc21[1], par$cutinc22[1],
        par$cutinc23[1])), Inf)
      
      eta1 <- par$mu1
      eta2 <- par$mu2
      
      rho <- par$rho
      
      d <- loglik(eta1 = eta1, eta2 = eta2, rho = rho, alpha1 = alpha1,
        alpha2 = alpha2, y = y, log = log)
      
      return(d)
    }
  )
#  f$initialize <- list(
#    "rho" = function(y, ...) rep(0.5, length(y))
#  )
  class(f) <- "gamlss2.family"
  return(f)
}

