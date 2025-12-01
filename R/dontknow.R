## Helper functions.
inc2cut <- function(x)
{
  if(length(x) > 1L) {
    x[-1L] <- x[1L] + cumsum(exp(x[-1L]))
  }
  return(x)
}

## Bivariate standard normal pdf at (h, k) with correlation rho.
DK_bvn_pdf <- function(h, k, rho) {
  if(!is.finite(h) || !is.finite(k)) return(0)

  r <- max(min(rho, 0.999999999), -0.999999999)
  omr2 <- 1 - r^2
  if(omr2 <= 0) return(0)

  quad <- (h^2 - 2 * r * h * k + k^2) / omr2
  pref <- 1 / (2 * pi * sqrt(omr2))
  pref * exp(-0.5 * quad)
}

## dell phi2 / dell rho at (h, k; rho).
DK_bvn_pdf_drho <- function(h, k, rho) {
  if(!is.finite(h) || !is.finite(k)) return(0)

  r <- max(min(rho, 0.999999999), -0.999999999)
  omr2 <- 1 - r^2
  if(omr2 <= 0) return(0)

  cval <- h^2 - 2 * r * h * k + k^2

  ## d/dρ log φ₂(h,k;ρ)
  dlogphi <- r / omr2 +
    (h * k * omr2 - r * cval) / (omr2^2)

  phi <- DK_bvn_pdf(h, k, r)
  phi * dlogphi
}

## rhogit inverse-link derivatives expressed in terms of rho
DK_rhogit_d1 <- function(rho) {
  r <- pmax(pmin(rho, 0.999999999), -0.999999999)
  omr2 <- pmax(1 - r^2, .Machine$double.eps)
  omr2^(3/2)          # d rho / d eta
}

DK_rhogit_d2 <- function(rho) {
  r <- pmax(pmin(rho, 0.999999999), -0.999999999)
  omr2 <- pmax(1 - r^2, .Machine$double.eps)
  -3 * r * omr2^2     # d^2 rho / d eta^2
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

  ## Derivatives for rho via C (score wrt eta_rho, Fisher hess wrt eta_rho)
  f$score <- list()
  f$hess  <- list()

  f$score$rho <- function(y, par, ...) {
    y <- as.matrix(y)

    ## alpha matrix with increasing cuts (same as in pdf/logLik)
    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
    }

    .Call(
      "score_rho_dontknow",
      as.numeric(par$mu1),    # Eta1 on mean scale
      as.numeric(par$mu2),    # Eta2 on mean scale
      as.numeric(par$rho),    # rho on parameter scale
      alpha,                  # n x mA alpha matrix (already monotone)
      y,                      # n x 2 response (y1, y2) with y2 in {0,...}
      PACKAGE = "dontknow"
    )
  }

  f$hess$rho <- function(y, par, ...) {
    ## y not needed for Fisher info (we integrate over categories), but
    ## signature is fixed, so we still accept it.

    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
    }

    .Call(
      "hess_rho_dontknow",
      as.numeric(par$mu1),    # Eta1
      as.numeric(par$mu2),    # Eta2
      as.numeric(par$rho),    # rho (parameter scale)
      alpha,                  # n x mA alpha
      PACKAGE = "dontknow"
    )
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

sim_DK <- function(n = 1000, rho = 0.5)
{
  stopifnot(requireNamespace("mvtnorm"))

  d <- data.frame(
    "x1" = runif(n, -3, 3),
    "x2" = runif(n, -3, 3),
    "x3" = runif(n, -3, 3)
  )

  mu1 <- sin(d$x1) * 2
  mu2 <- d$x2^2 - 3
  rl <- make.link2("rhogit")
  if(is.null(rho)) {
    rho <- rl$linkinv(cos(d$x3))
  } else {
    rho <- rep(rho, length.out = n)
  }

  alpha1 <- c(-Inf, 0, Inf)
  alpha2 <- c(-Inf, -1, 0, 1, Inf)

  yhelp <- NULL
  for(i in 1:n) {
    Sigma <- matrix(c(1, rho[i], rho[i], 1), 2, 2)
    yhelp <- rbind(yhelp, mvtnorm::rmvnorm(1, mean = c(0, 0), sigma = Sigma))
  }

  y1star <- mu1 + yhelp[, 1]
  y2star <- mu2 + yhelp[, 2]

  d$y1 <- cut(y1star, alpha1, labels = FALSE) - 1
  d$y2 <- cut(y2star, alpha2, labels = FALSE) - 1

  ## new response
  d$Y <- cbind(d$y1, d$y2)

  return(d)
}

