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

  ## convenience to build monotonized alpha matrix
  build_alpha <- function(par) {
    a <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      a[, -1L] <- t(apply(a[, -1L], 1L, inc2cut))
    }
    a
  }

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
    probs <- f$probabilities(par)          ## n x (K_cat+1)
    n     <- nrow(probs)
    K_cat <- ncol(probs) - 1L              ## Y2 categories 0, ..., K_cat-1

    ## coerce y
    if(length(y) == 1L)
      y <- rep.int(y, n)
    y <- as.integer(y)

    if(any(y < 0L | y > (K_cat - 1L)))
      stop("y must be in 0 ...", K_cat - 1L)

    ## cumulative along order: DK, Y2=0, Y2=1, ...
    cprobs <- t(apply(probs, 1L, cumsum))

    ## index: DK = col 1, Y2 = c -> col c+2
    ans <- cprobs[cbind(seq_len(n), y + 1L)]

    if(!lower.tail)
      ans <- 1 - ans + probs[cbind(seq_len(n), y + 1L)]
    if(log.p)
      ans <- log(ans)

    ans
  }

  f$probabilities <- function(par) {
    eta1  <- par$mu1
    eta2  <- par$mu2
    rho   <- par$rho
    alpha <- build_alpha(par)   ## ordered thresholds: alpha1 (DK), alpha2, ..., alpha_k (Y2 cuts)

    stopifnot(requireNamespace("mvtnorm"))

    n     <- length(eta1)
    Tcuts <- ncol(alpha) - 1L          ## number of finite cuts for Y2
    K_cat <- Tcuts + 1L                ## number of Y2 categories: 0, ..., K_cat-1

    ## matrix: col1 = DK, cols 2..(K_cat+1) = Y2 = 0, ..., K_cat-1
    out <- matrix(NA_real_, n, K_cat + 1L)

    for(i in seq_len(n)) {
      r <- rho[if(length(rho) == 1L) 1L else i]
      r <- max(min(r, 0.999999999), -0.999999999)

      ## DK probability: P(Y1 = 1) = 1 - Phi(alpha1 - eta1)
      Aup <- alpha[i, 1L] - eta1[i]
      out[i, 1L] <- pnorm(Aup, lower.tail = FALSE)

      Sigma <- matrix(c(1, r, r, 1), 2, 2)

      ## finite cuts for Y2 on latent scale (Z2), for this i
      ## alpha2..alpha_k, so length(Tcuts)
      cuts <- alpha[i, -1L] - eta2[i]

      ## non-DK categories for Y2 = 0, ..., K_cat-1
      for(c in 0:(K_cat - 1L)) {
        B1 <- if(c == 0L)          -Inf else cuts[c]       ## lower bound
        B2 <- if(c == (K_cat - 1L)) Inf else cuts[c + 1L]  ## upper bound

        p_up <- if(is.infinite(B2) && B2 > 0) {
          ## P(Z1 <= Aup) when Z2 <= +Inf
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

    ## guard against small numerical drift
    rs <- rowSums(out)
    ok <- rs > 0 & is.finite(rs)
    out[ok, ] <- out[ok, ] / rs[ok]

    ## column names: "DK", "Y2.0", "Y2.1", ..., "Y2.(K_cat-1)"
    colnames(out) <- c("DK", paste0("Y2.", seq_len(ncol(out) - 1L) - 1L))

    out
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

  ## Derivatives
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

  ## score wrt eta_mu1
  f$score$mu1 <- function(y, par, ...) {
    y <- as.matrix(y)

    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
    }

    .Call(
      "score_mu1_dontknow",
      as.numeric(par$mu1),    # Eta1
      as.numeric(par$mu2),    # Eta2
      as.numeric(par$rho),    # rho
      alpha,                  # n x mA alpha
      y,                      # n x 2 Y (y1,y2), y2 starts at 0
      PACKAGE = "dontknow"
    )
  }

  ## score wrt eta_mu2
  f$score$mu2 <- function(y, par, ...) {
    y <- as.matrix(y)

    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
    }

    .Call(
      "score_mu2_dontknow",
      as.numeric(par$mu1),
      as.numeric(par$mu2),
      as.numeric(par$rho),
      alpha,
      y,
      PACKAGE = "dontknow"
    )
  }

  ## Hessian (Fisher info) wrt eta_mu1
  f$hess$mu1 <- function(y, par, ...) {
    ## y is not needed for Fisher, covariates/parameters are enough.

    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
    }

    .Call(
      "hess_mu1_dontknow",
      as.numeric(par$mu1),
      as.numeric(par$mu2),
      as.numeric(par$rho),
      alpha,
      PACKAGE = "dontknow"
    )
  }

  ## Hessian (Fisher info) wrt eta_mu2
  f$hess$mu2 <- function(y, par, ...) {

    alpha <- do.call("cbind", par[alpha_names])
    if(k > 2L) {
      alpha[, -1L] <- t(apply(alpha[, -1L], 1L, inc2cut))
    }

    .Call(
      "hess_mu2_dontknow",
      as.numeric(par$mu1),
      as.numeric(par$mu2),
      as.numeric(par$rho),
      alpha,
      PACKAGE = "dontknow"
    )
  }

  for(jj in seq_along(alpha_names)) {
    par_name <- alpha_names[jj]       # "alpha1", "alpha2", ...
    j_idx <- jj                       # 1-based column index

    ## score wrt alpha_j
    f$score[[par_name]] <- local({
      j <- j_idx
      function(y, par, ...) {
        y <- as.matrix(y)
        alpha <- build_alpha(par)
        .Call(
          "score_alpha_dontknow",
          as.numeric(par$mu1),
          as.numeric(par$mu2),
          as.numeric(par$rho),
          alpha,
          y,
          as.integer(j),
          PACKAGE = "dontknow"
        )
      }
    })

    ## hessian (Fisher) wrt alpha_j
    f$hess[[par_name]] <- local({
      j <- j_idx
      function(y, par, ...) {
        alpha <- build_alpha(par)
        .Call(
          "hess_alpha_dontknow",
          as.numeric(par$mu1),
          as.numeric(par$mu2),
          as.numeric(par$rho),
          alpha,
          as.integer(j),
          PACKAGE = "dontknow"
        )
      }
    })
  }

  f$alpha <- build_alpha

  f$residuals <- function(object, DK = TRUE, ...) {
    if(DK) {
      return(rqres_DK_indicator(object, ...))
    } else {
      return(rqres_DK_ordinal(object, ...))
    }
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

rqres_DK_indicator <- function(object, ...)
{
  fam <- family(object)
  if(!identical(fam$family, "DK"))
    stop("DK family required.")

  mf <- model.frame(object)
  Y  <- mf$Y
  par   <- predict(object)
  probs <- fam$probabilities(par)

  y1 <- Y[, 1L]
  p1 <- probs[, 1L] ## P(DK)

  ## probabilities for Y1=0 and Y1=1
  p0 <- 1 - p1

  F_upper <- ifelse(y1 == 0L, p0, 1)
  F_lower <- ifelse(y1 == 0L, 0, p0)

  u <- runif(length(y1), F_lower, F_upper)
  qnorm(u)
}

rqres_DK_ordinal <- function(object, ...)
{
  fam <- family(object)
  if(!identical(fam$family, "DK"))
    stop("DK family required.")

  mf   <- model.frame(object)
  Y    <- mf$Y
  par  <- predict(object)
  probs <- fam$probabilities(par)

  y1 <- Y[, 1L]
  y2 <- Y[, 2L]

  ## only non-DK observations
  keep <- (y1 == 0L)
  if(!any(keep))
    stop("No non-DK observations.")

  probs_ndk <- probs[keep, , drop = FALSE]
  y2_ndk    <- y2[keep]

  ## probs_ndk: col1 = DK, cols 2..(K+1) = Y2 = 0..K-1
  K <- ncol(probs_ndk) - 1L

  ## check that y2 is in the supported range 0, ..., K-1
  if(any(y2_ndk < 0L | y2_ndk > (K - 1L))) {
    bad_vals <- sort(unique(y2_ndk[y2_ndk < 0L | y2_ndk > (K - 1L)]))
    stop(
      "y2 has values outside the supported range 0..", K - 1L,
      " implied by DK(k).\n",
      "Offending values: ", paste(bad_vals, collapse = ", "), "\n",
      "Either recode y2 to 0..", K - 1L,
      " or adjust the DK(k) specification so that all categories are represented."
    )
  }

  ## conditional probabilities for Y2 = c given Y1 = 0
  p_dk  <- probs_ndk[, 1L]
  p_y2  <- probs_ndk[, -1L, drop = FALSE] ## cols for Y2 = 0, ..., K-1
  p_tot <- 1 - p_dk                       ## = P(Y1 = 0)

  ## we can only define conditional probs where P(Y1=0) > 0
  ok_row <- is.finite(p_tot) & (p_tot > 0)

  if(!any(ok_row)) {
    warning("No non-DK observations with positive predicted P(Y1=0); returning NA.")
    return(rep(NA_real_, length(y1)))
  }

  ## conditional probs only for "ok" rows
  p_cond <- matrix(NA_real_, nrow(p_y2), ncol(p_y2))
  p_cond[ok_row, ] <- p_y2[ok_row, , drop = FALSE] / p_tot[ok_row]

  ## cumulative over Y2 categories 0..K-1
  cprobs_cond <- t(apply(p_cond, 1L, cumsum))

  z <- y2_ndk

  ## indices of rows where we have a proper conditional distribution
  ok_z <- ok_row & !is.na(z)

  ## F_upper for all ok_z
  F_upper_ndk <- rep(NA_real_, length(z))
  idx_up <- which(ok_z)
  if(length(idx_up) > 0L) {
    F_upper_ndk[idx_up] <- cprobs_cond[cbind(idx_up, z[idx_up] + 1L)]
  }

  ## F_lower: do NOT use ifelse with matrix indexing
  F_lower_ndk <- rep(NA_real_, length(z))
  ## z == 0: lower bound is 0
  idx0 <- which(ok_z & z == 0L)
  if(length(idx0) > 0L)
    F_lower_ndk[idx0] <- 0

  ## z > 0: lower bound is cumulative up to category (z-1)
  idx_pos <- which(ok_z & z > 0L)
  if(length(idx_pos) > 0L)
    F_lower_ndk[idx_pos] <- cprobs_cond[cbind(idx_pos, z[idx_pos])]

  ## clamp to [0, 1] and fix potential tiny numerical inversions
  F_upper_ndk <- pmin(pmax(F_upper_ndk, 0), 1)
  F_lower_ndk <- pmin(pmax(F_lower_ndk, 0), 1)

  bad_bounds <- !is.na(F_lower_ndk) & !is.na(F_upper_ndk) &
    (F_lower_ndk > F_upper_ndk)
  if(any(bad_bounds)) {
    tmp <- F_lower_ndk[bad_bounds]
    F_lower_ndk[bad_bounds] <- F_upper_ndk[bad_bounds]
    F_upper_ndk[bad_bounds] <- tmp
  }

  ## simulate u only where we have valid bounds
  good <- !is.na(F_lower_ndk) & !is.na(F_upper_ndk)
  u_ndk <- rep(NA_real_, length(z))
  if(any(good))
    u_ndk[good] <- runif(sum(good), F_lower_ndk[good], F_upper_ndk[good])

  r_ndk <- ifelse(is.na(u_ndk), NA_real_, qnorm(u_ndk))

  ## plug back into full residual vector (non-DK only)
  res <- rep(NA_real_, length(y1))
  res[keep] <- r_ndk

  res
}

