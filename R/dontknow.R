## Helper functions.
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

      ## Second dimension bounds depend on observed y2 in {0,1,2,3,...}.
      B2 <- alpha2[y2[i] + 2L] - eta2[i] ## upper: alpha_{2, c} - eta2  (c = y2[i]+1)
      B1 <- alpha2[y2[i] + 1L] - eta2[i] ## lower: alpha_{2, c-1} - eta2

      ## Prob(A, B2; rho) - Prob(A, B1; rho)
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
  ## Parameter names in a guaranteed order.
  alpha_names <- c("alpha1", paste0("alpha", 2:k))
  names <- c("mu1", "mu2", "rho", alpha_names)
  links <- c("identity", "identity", "rhogit", rep("identity", length(alpha_names)))

  f <- list(
    "family" = "DK",
    ## Parameterization: alpha1 is the binary cut; alpha2 ... alphak are the ordinal cuts.
    "names"  = names,
    "links"  = links,
    "pdf" = function(y, par, log = FALSE) {
      ## ensure y2 in {0,1,...,K-1}
      y2min <- min(y[, 2L], na.rm = TRUE)
      if(y2min != 0L) {
        y[, 2L] <- y[, 2L] - y2min
      }

      ## Build the ordinal cut matrix and ensure ordering.
      alpha <- do.call("cbind", par[alpha_names])
      if(k > 2L) {
        alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
      }

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

  f$cdf <- function(y, par, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- length(par$mu1)

    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L)
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))

    probs <- cdf_dontknow_R(par$mu1, par$mu2, par$rho, alpha)

    ## probs: n x (K+1), col1=dk, col(c+2)=category c
    ## build cumulative along that order:
    cprobs <- t(apply(probs, 1L, cumsum))  ## n x (K+1)

    ## coerce
    if(length(y) == 1L)
      y <- rep.int(y, n)
    y <- as.integer(y)
    if(any(y < 0L | y > (k-1L))) stop("y must be in 0..", k-1L)

    ans <- cprobs[cbind(seq_len(n), y + 1L)]

    if(!lower.tail)
      ans <- 1 - ans + probs[cbind(seq_len(n), y + 1L)]
    if(log.p)
      ans <- log(ans)

    ans
  }

  f$probabilities <- function(par) {
    cdf_dontknow(par)
  }

  f$valid.response <- function(x) {
    if(is.factor(x) | is.character(x)) 
      stop("the response should be integer/numeric!")
    if(!all(range(x[, 1L]) %in% c(0, 1)))
      stop("response y1 must be 0/1!")
    if(min(x[, 2L]) != 0)
      stop("response y2 must be 0, 1, 2, 3, ...!")
    return(TRUE)
  }

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

pbvnorm_miwa <- function(h, k, rho, steps = 128L) {
  .Call("pbvnorm_miwa", as.numeric(h), as.numeric(k), as.numeric(rho), as.integer(steps))
}

## CDF function.
cdf_dontknow <- function(par)
{
  eta1 <- par$mu1
  eta2 <- par$mu2
  rho <- par$rho
  alpha <- as.matrix(par[, grep("alpha", names(par))])

  stopifnot(requireNamespace("mvtnorm"))
  n <- length(eta1)
  K <- ncol(alpha) - 1L
  out <- matrix(NA_real_, n, K + 1L)

  for(i in seq_len(n)) {
    r <- rho[if(length(rho) == 1L) 1L else i]
    r <- max(min(r,  0.999999999), -0.999999999)

    Aup <- alpha[i, 1L] - eta1[i]
    out[i, 1L] <- pnorm(Aup, lower.tail = FALSE)

    for(c in 0:(K-1L)) {
      B1 <- if(c == 0L) -Inf else alpha[i, c+1L] - eta2[i]
      B2 <- if(c+1L >= K)  Inf else alpha[i, c+2L] - eta2[i]

      Sigma <- matrix(c(1, r, r, 1), 2, 2)

      p_up <- if(is.infinite(B2) && B2 > 0) {
        pnorm(Aup)
      } else {
        as.numeric(mvtnorm::pmvnorm(
          lower = c(-Inf, -Inf),
          upper = c(Aup, B2),
          mean  = c(0, 0),
          sigma = Sigma,
          keepAttr = FALSE
        ))
      }
      p_lo <- if(is.infinite(B1) && B1 < 0) {
        0
      } else {
        as.numeric(mvtnorm::pmvnorm(
          lower = c(-Inf, -Inf),
          upper = c(Aup, B1),
          mean  = c(0, 0),
          sigma = Sigma,
          keepAttr = FALSE
        ))
      }
      pc <- p_up - p_lo
      out[i, c + 2L] <- if(pc > 0) pc else 0
    }
  }
  out
}

