\name{DK}
\alias{DK}
\alias{dontknow}
\alias{sim_DK}

\title{Don't Know (DK) Family}

\description{
  Defines the don't know family.
}

\usage{
## family function
DK(k = 4, useC = TRUE)

## simulation function for DK(4)
sim_DK(n = 1000, rho = NULL, shift = 0)
}

\arguments{
  \item{k}{Integer, the number of levels of \code{y2}. Default is \code{k = 4}.}
  \item{useC}{Logical, should the C version of the log-likelihood be used?}
  \item{n}{Integer, number of observations to simulate.}
  \item{rho}{Numeric, the correlation parameter. If \code{rho = NULL}
    \code{rho} is a function of covariate \code{x3}.}
  \item{shift}{Numeric shift added to the latent mean predictor \code{mu1}
    in the simulation design.}
}

\details{
  \code{DK()} defines a \code{"gamlss2.family"} object for models with an
  explicit \dQuote{don't know} response option. The response is represented
  as a two-column matrix \code{Y = cbind(y1, y2)} where \code{y1} is a
  binary indicator for choosing the \dQuote{don't know} option and
  \code{y2} is an ordinal response for the non-\dQuote{don't know}
  categories, coded as \code{0, 1, \dots, k - 1}.

  \code{sim_DK()} generates synthetic data from the same latent-variable
  construction used in the examples, returning covariates, latent effects,
  latent responses, and the bivariate response matrix \code{Y}.
}

\value{
  \code{DK()} returns a \code{"gamlss2.family"} object to be used with
  \code{\link{gamlss2}}.

  \code{sim_DK()} returns a data frame containing simulated covariates,
  latent predictors, latent responses, observed responses \code{y1} and
  \code{y2}, and the combined bivariate response \code{Y}.
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
alpha <- family(b)$alpha(par)
head(alpha)

## Plot ystars.
probs <- family(b)$probabilities(par)
col_fit <- apply(probs, 1, which.max)
col_truth <- d$y2 + 1L
col_truth[d$y1 == 1L] <- 5L

par(mfrow = c(1, 2))
plot(y2star ~ y1star, data = d, col = col_fit, main = "Estimate")
abline(h = alpha, lty = 2, col = "lightgray")
plot(y2star ~ y1star, data = d, col = col_truth, main = "Truth")

## DK calibration plot.
calibration(b)

## randomized quantile residuals
r_DK <- residuals(b, DK = TRUE)
r_cat <- residuals(b, DK = FALSE)

par(mfrow = c(2, 3), oma = c(0, 0, 0, 4))
plot(r_DK, spar = FALSE)
mtext("DK residuals", side = 4, line = 2, font = 2)
plot(r_cat, spar = FALSE)
mtext("CAT residuals", side = 4, line = 2, font = 2)

## simulation study
set.seed(1328)
d <- sim_DK(2000)

## formula
f <- Y ~ -1 + s(x1) + la(id1,type=4) + la(id4,type=2) + la(id5,type=2) |
  -1 + s(x2) + la(id2,type=4) + la(id6,type=2) + la(id7,type=2) |
   1 + s(x3) + la(id3,type=4) + la(id8,type=2) + la(id9,type=2)

## estimate model
b <- gamlss2(f, family = DK(4), data = d)

## plot effects
plot(b)

## coefficient paths
plot_lasso(b, which = "coefficients", term = "mu1")
plot_lasso(b, which = "coefficients", term = "mu2")
plot_lasso(b, which = "coefficients", term = "rho")
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}
