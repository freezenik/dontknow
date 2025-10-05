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
  y1 <- y[, 1]
  y2 <- y[, 2]
  
  mu <- c(0, 0)
  
  ll <- rep(0, length(y1))
  
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2)
  
  index.ones <- which(y1 == 1)
  index.zeros <- which(y1 == 0)
  
  ll[index.ones] <- log(1 - pnorm(alpha1[2] - eta1[index.ones]))
  
  Sigma[1, 2] <- Sigma[2, 1] <- rho[1]
  chSigma <- chol(Sigma)
  chSigma <- as.ltMatrices(t(chSigma))

  lower <- rbind(rep(-Inf, length(index.zeros)),
    alpha2[y2[index.zeros]+1] - eta2[index.zeros])
  upper <- rbind(alpha1[2] - eta1[index.zeros],
    alpha2[y2[index.zeros]+2] - eta2[index.zeros])

  ll[index.zeros] <- lpmvnorm(
    lower = lower,
    upper = upper,
    mean = mu,
    chol = chSigma,
    logLik = FALSE,
    M = 250
  )
  
  if(!log)
    ll <- exp(ll)
  
  return(ll)
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
      
      rho <- invrhogit(par$rho)
      
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
      
      rho <- invrhogit(par$rho)
      
      d <- loglik(eta = eta1, eta2 = eta2, rho = rho, alpha1 = alpha1,
        alpha2 = alpha2, y = y, log = log)
      
      return(d)
    }
  )
  class(f) <- "gamlss2.family"
  return(f)
}

