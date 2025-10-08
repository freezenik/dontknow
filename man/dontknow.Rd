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
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}, \code{\link[MASS]{polr}}
}

\examples{\dontrun{## simulated data
library("gamlss2")
library("mvtnorm")

set.seed(1328)
rho <- 0.5

## number of observations
n <- 500

d <- data.frame(
  "x1" = runif(n, -3, 3),
  "x2" = runif(n, -3, 3)
)

mu1 <- sin(d$x1)
mu2 <- 0.5 * d$x2

alpha1 <- c(-Inf, 0, Inf)
alpha2 <- c(-Inf, -1, 0, 1, Inf)

Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
yhelp <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)

y1star <- mu1 + yhelp[, 1]
y2star <- mu2 + yhelp[, 2]

d$y1 <- cut(y1star, alpha1, labels = FALSE) - 1
d$y2 <- cut(y2star, alpha2, labels = FALSE) - 1

## new response
d$Y <- cbind(d$y1, d$y2)

## important, remove intercepts!
f <- Y ~ -1 + s(x1) | -1 + s(x2)

## estimate model
b <- gamlss2(f, family = DK, data = d)

## plot estimated effects
plot(b)

## extract estimated parameters
par <- predict(b, type = "parameter")

## check rho
head(par$rho)

## compare mu1 and mu2.
par(mfrow = c(1, 2))

plot(d$x1, y1star, pch = 16,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.1))
i <- order(d$x1)
lines(mu1[i] ~ d$x1[i], col = 2, lwd = 2)
lines(par$mu1[i] ~ d$x1[i], col = 4, lwd = 2)

plot(d$x2, y2star, pch = 16,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.1))
i <- order(d$x2)
lines(mu2[i] ~ d$x2[i], col = 2, lwd = 2)
lines(par$mu2[i] ~ d$x2[i], col = 4, lwd = 2)
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

