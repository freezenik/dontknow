\name{dk4}
\alias{DK}
\alias{DK}
\alias{dontknow}

\title{Don't Know (DK) Family}

\description{
  Defines the don't know family.
}

\usage{
DK(k = 4, useC = TRUE)
}

\arguments{
  \item{k}{Integer, the number of levels of y2. Defaul is \code{k = 4}.}
  \item{useC}{Logical, should the C version of the log-likelihood be used?}
}

\details{
  Details comming soon!
}

\value{
  A \code{"gamlss2.family"} object to be used with \code{\link{gamlss2}}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}
}

\examples{\dontrun{## simulated data
library("gamlss2")
library("mvtnorm")

set.seed(1328)

## number of observations
n <- 1000

d <- data.frame(
  "x1" = runif(n, -3, 3),
  "x2" = runif(n, -3, 3),
  "x3" = runif(n, -3, 3)
)

mu1 <- sin(d$x1) * 2
mu2 <- d$x2^2 - 3
rl <- make.link2("rhogit")
rho <- rl$linkinv(cos(d$x3))

alpha1 <- c(-Inf, 0, Inf)
alpha2 <- c(-Inf, -1, 0, 1, Inf)

yhelp <- NULL
for(i in 1:n) {
  Sigma <- matrix(c(1, rho[i], rho[i], 1), 2, 2)
  yhelp <- rbind(yhelp, rmvnorm(1, mean = c(0, 0), sigma = Sigma))
}

y1star <- mu1 + yhelp[, 1]
y2star <- mu2 + yhelp[, 2]

d$y1 <- cut(y1star, alpha1, labels = FALSE) - 1
d$y2 <- cut(y2star, alpha2, labels = FALSE) - 1

## new response
d$Y <- cbind(d$y1, d$y2)

## important, remove intercepts!
f <- Y ~ -1 + s(x1) | -1 + s(x2) | s(x3)

## estimate model, y2 with 4 levels
b <- gamlss2(f, family = DK(4), data = d)

## plot estimated effects
plot(b)

## extract estimated parameters
par <- predict(b, type = "parameter")

## compare mu1, mu2 and rho
par(mfrow = c(1, 3))

plot(d$x1, y1star, pch = 16,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.1),
  main = "mu1")
i <- order(d$x1)
lines(mu1[i] ~ d$x1[i], col = 2, lwd = 2)
lines(par$mu1[i] ~ d$x1[i], col = 4, lwd = 2)

plot(d$x2, y2star, pch = 16,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.1),
  main = "mu2")
i <- order(d$x2)
lines(mu2[i] ~ d$x2[i], col = 2, lwd = 2)
lines(par$mu2[i] ~ d$x2[i], col = 4, lwd = 2)

i <- order(d$x3)
plot(d$x3[i], par$rho[i], type = "l", col = 2, lwd = 2,
  ylim = range(rho, par$rho), main = "rho")
lines(rho[i] ~ d$x3[i], col = 4, lwd = 2)

## check cuts
inc2cut <- function(x)
{
  if(length(x) > 1L) {
    x[-1L] <- x[1L] + cumsum(exp(x[-1L]))
  }
  return(x)
}

alpha <- t(apply(par[, 5:7], 1, inc2cut))
head(alpha)
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

