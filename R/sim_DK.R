sim_DK <- function(n = 1000, rho = NULL, shift = 0)
{
  stopifnot(requireNamespace("mvtnorm"))

  d <- data.frame(
    "x1" = runif(n, -3, 3),
    "x2" = runif(n, -3, 3),
    "x3" = runif(n, -3, 3)
  )

  for(j in 1:2) {
    d[[paste0("id", j)]] <- factor(sample(1:5, size = n, replace = TRUE))
  }
  d[[paste0("id", 3)]] <- factor(sample(1:7, size = n, replace = TRUE))
  for(j in c(4, 6, 8)) {
    d[[paste0("id", j)]] <- factor(sample(1:4, size = n, replace = TRUE))
  }
  for(j in c(5, 7, 9)) {
    d[[paste0("id", j)]] <- factor(sample(1:6, size = n, replace = TRUE))
  }

  b1 <- c(0, 1.5, -1.2, -0.5, 1.5)
  b2 <- c(0, 1.2, 1.2, -2, -1.5)
  b3 <- c(0, 0.5, -1, -2, -2, 1, 0.5)

  X1 <- model.matrix(~ id1, data = d)
  X2 <- model.matrix(~ id2, data = d)
  X3 <- model.matrix(~ id3, data = d)

  d$fx1 <- sin(d$x1) * 2
  d$fx2 <- d$x2^2 / 4.5
  d$fx3 <- if(is.null(rho)) d$x3 * 2/3 else rep(0, n)
  d$fx1 <- d$fx1 - mean(d$fx1)
  d$fx2 <- d$fx2 - mean(d$fx2)
  d$fx3 <- d$fx3 - mean(d$fx3)

  d$fid1 <- drop(X1 %*% b1)
  d$fid2 <- drop(X2 %*% b2)
  d$fid3 <- if(is.null(rho)) drop(X3 %*% b3) else rep(0, n)

  d$mu1 <- d$fx1 + d$fid1 + shift
  d$mu2 <- d$fx2 + d$fid2

  if(is.null(rho)) {
    d$rho <- d$fx3 + d$fid3
  } else {
    d$rho <- rep(rho, length.out = n)
  }

  rl <- make.link2("rhogit")
  d$rho <- rl$linkinv(d$rho)

  alpha1 <- c(-Inf, 0, Inf)
  alpha2 <- c(-Inf, -1, 0, 1, Inf)

  yhelp <- NULL
  for(i in 1:n) {
    Sigma <- matrix(c(1, d$rho[i], d$rho[i], 1), 2, 2)
    yhelp <- rbind(yhelp, mvtnorm::rmvnorm(1, mean = c(0, 0), sigma = Sigma))
  }

  d$y1star <- d$mu1 + yhelp[, 1]
  d$y2star <- d$mu2 + yhelp[, 2]

  d$y1 <- cut(d$y1star, alpha1, labels = FALSE) - 1
  d$y2 <- cut(d$y2star, alpha2, labels = FALSE) - 1

  d$Y <- cbind(d$y1, d$y2)

  return(d)
}

