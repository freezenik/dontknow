## Helper functions.
invrhogit <- function(eta)
{
  return(eta / sqrt(1 + eta^2))
}

rhogit <- function(mu)
{
  return(mu / sqrt(1 - mu^2))
}

cut2inc <- function(x)
{
  if(length(x) > 1) {
    x[-1L] <- log(diff(x))
  }
  return(x)  
}

inc2cut <- function(x)
{
  if(length(x) > 1L) {
    x[-1L] <- x[1L] + cumsum(exp(x[-1L]))
  }
  return(x)
}

## logLik function.
logLik_dontknow <- function(eta1, eta2, rho, alpha, y, log = TRUE)
{
  stopifnot(requireNamespace("mvtnorm"))

  y1 <- y[, 1L]
  y2 <- y[, 2L]
  n  <- length(y1)

  ll <- numeric(n)
  rho <- rep(rho, length.out = n)

  for(i in seq_len(n)) {
    ## Correlation matrix for this obs.
    Sigma <- matrix(c(1, rho[i], rho[i], 1), 2, 2)

    if(y1[i] == 1L) {
      ## P(Y1 = 1) = 1 - Phi(alpha1 - eta1).
      if(log) {
        ll[i] <- pnorm(alpha[i, 1] - eta1[i], lower.tail = FALSE, log.p = TRUE)
      } else {
        ll[i] <- pnorm(alpha[i, 1] - eta1[i], lower.tail = FALSE)
      }
    } else {
      alpha2 <- c(-Inf, alpha[i, -1], Inf)

      ## y1 == 0: rectangle on (Z1, Z2).
      ## First dimension upper bound is always alpha1[i] - eta1[i].
      A  <- alpha[i, 1] - eta1[i]

      ## Second dimension bounds depend on observed y2 in {0,1,2,3}.
      B2 <- alpha2[y2[i] + 2L] - eta2[i] ## upper: alpha_{2, c} - eta2  (c = y2[i]+1)
      B1 <- alpha2[y2[i] + 1L] - eta2[i] ## lower: alpha_{2, c-1} - eta2

      ## Φ2(A, B2; ρ) - Φ2(A, B1; ρ)
      p_up <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = c(A, B2),
                               mean = c(0, 0), sigma = Sigma, keepAttr = FALSE)
      p_lo <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = c(A, B1),
                               mean = c(0, 0), sigma = Sigma, keepAttr = FALSE)
      p <- as.numeric(p_up) - as.numeric(p_lo)
      if(log) {
        ## Guard tiny numerical negatives.
        ll[i] <- if(p > 0) log(p) else -Inf
      } else {
        ll[i] <- if(p > 0) p else 0
      }
    }
  }

  return(ll)
}

## Don't-know family: explicit thresholds as parameters.
DK <- function(k = 4, useC = TRUE)
{
  names <- c("mu1", "mu2", "rho", "alpha1")
  names <- c(names, paste0("alpha", 2:k))
  links <- c("identity", "identity", "rhogit", rep("identity", k))

  f <- list(
    family = "DK",
    ## Parameterization: alpha1 is the binary cut; alpha2 ... alphak are the ordinal cuts.
    names  = names,
    links  = links,
    d = function(y, par, log = FALSE) {
      ## Build the ordinal cut matrix and ensure ordering.
      alpha <- do.call("cbind", par[grep("alpha", names(par))])
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))

     ll <- if(isTRUE(useC)) {
       logLik_dontknow_C
     } else {
       logLik_dontknow
     }

      ## Evaluate log-likelihood.
      d <- ll(
        eta1   = par$mu1,
        eta2   = par$mu2,
        rho    = par$rho,
        alpha  = alpha,
        y      = y,
        log    = log
      )

      return(d)
    }
  )
  class(f) <- "gamlss2.family"
  f
}

## Fast C version.
logLik_dontknow_C <- function(eta1, eta2, rho, alpha, y, log = TRUE) {
  .Call("logLik_dontknow",
        as.numeric(eta1),
        as.numeric(eta2),
        as.numeric(rho),
        as.matrix(alpha),
        as.matrix(y),
        as.logical(log), package = "dontknow")
}

