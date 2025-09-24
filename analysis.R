## Required packages.
library("readxl")
if(!("gamlss2" %in% installed.packages())) {
  install.packages("gamlss2",
    repos = c("https://gamlss-dev.R-universe.dev",
              "https://cloud.R-project.org"))
}
library("gamlss2")

## Functions.
source("functions.R")

## Create figures directory.
if(!file.exists("figures")) {
  dir.create("figures")
}

## Data.
d <- read_excel("data/Dati_bdI_DK.xls")
d <- as.data.frame(d)
names(d) <- gsub('\"', '', names(d), fixed = TRUE)

d <- with(d, data.frame(
  "age" = ETA_12, 
  "gender" = SEX,
  "marital" = STACIV_12,
  "geo" = AREA5_12,
  "status" = ASPRED_12,
  "hhs" = NCOMP_12,
  "y1" = VARRED_12DK,
  "y2" = VARRED_12
))

## After discussion with Niki, it is important that
## all factor variables are encoded as such.
d$gender <- d$gender - 1
d$gender <- as.factor(d$gender)
d$marital <- as.factor(d$marital)
d$geo <- as.factor(d$geo)
d$status <- as.factor(d$status)
d$hhs <- as.factor(d$hhs)

## Attention: in the current way Thomas and I have
## implemented the log-likelihood function,
## the ordinal response y has to be in {0, ..., C - 1}!!!
d$y2 <- d$y2 - 1

## Estimation in gamlss2.
## Response.
d$Y <- cbind(d$y1, d$y2)

## Formula.
## type = 1 -> simple lasso
## type = 2 -> group lasso
## type = 3 -> ordinal fused lasso
## type = 4 -> nominal fused lasso
f <- Y ~ -1 + s(age) + la(gender,type=2) + la(marital,type=2) +
  la(geo,type=4) + la(status,type=2) + la(hhs,type=3) | .

## Estimate model.
if(!file.exists("model.rds")) {
  b <- gamlss2(f, family = dk3, data = d, eps = .00001)
  saveRDS(b, file = "model.rds")
} else {
  b <- readRDS("model.rds")
}

## Figures.
pdf(file = "figures/model_paths.pdf", width = 14, height = 9)

ci <- 0.7
lwd <- 2

par(mfrow = c(2, 5), mar = c(4, 4, 4, 3), oma = c(0, 4, 0, 0))

plot_lasso(b, which = "coefficients",
  term = "mu1.la(gender", spar = FALSE,
  zoom = c(20, 1), main = "gender",
  info = TRUE, cex.info = ci, lwd = lwd)

mtext(bquote(eta[1]), side = 2, line = 5, font = 2, cex = 1.5)

plot_lasso(b, which = "coefficients",
  term = "mu1.la(marital", spar = FALSE,
  zoom = c(20, 5), main = "marital",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu1.la(geo", spar = FALSE,
  zoom = c(20, 1), main = "geo",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu1.la(status", spar = FALSE,
  zoom = c(9, 1), main = "status",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu1.la(hhs", spar = FALSE,
  zoom = c(7, 2), main = "hhs",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(gender", spar = FALSE,
  zoom = c(10, 1), main = "gender",
  info = TRUE, cex.info = ci, lwd = lwd)

mtext(bquote(eta[2]), side = 2, line = 5, font = 2, cex = 1.5)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(marital", spar = FALSE,
  zoom = c(20, 5), main = "marital",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(geo", spar = FALSE,
  zoom = c(9, 1), main = "geo",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(status", spar = FALSE,
  zoom = c(6, 6), main = "status",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(hhs", spar = FALSE,
  zoom = c(7, 2), main = "hhs",
  info = TRUE, cex.info = ci, lwd = lwd)

dev.off()

pdf(file = "figures/age_effects.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(b, model = "mu1", term = "age", spar = FALSE, ylab = "f(age)")
mtext(bquote(eta[1]), side = 3, line = 1, font = 2, cex = 1.5)

plot(b, model = "mu2", term = "age", spar = FALSE, ylab = "f(age)")
mtext(bquote(eta[2]), side = 3, line = 1, font = 2, cex = 1.5)

dev.off()

## Simulated data.
library("mvtnorm")

set.seed(1328)
rho <- 0.5

## Number of observations.
n <- 500

## Simulate features with x1 and x5 being signal,
## and the rest noise variables
x1 <- c(1:5)
x1 <- as.factor(sample(x1, n, replace = TRUE))
x2 <- c(1:4)
x2 <- as.factor(sample(x2, n, replace = TRUE))
x3 <- c(1:6)
x3 <- as.factor(sample(x3, n, replace = TRUE))
x4 <- c(1:5)
x4 <- as.factor(sample(x4, n, replace = TRUE))
x5 <- c(1:4)
x5 <- as.factor(sample(x5, n, replace = TRUE))
x6 <- c(1:6)
x6 <- as.factor(sample(x6, n, replace = TRUE))

X1 <- model.matrix(~ 1 + x1)
X4 <- model.matrix(~ 1 + x4)

beta_1 <- c(0, 0.5, -0.4, 0.3, -0.4)
beta_2 <- c(0, -0.4, 0.4, 0.5, -0.5)

mu1 <- X1 %*% beta_1
mu2 <- X4 %*% beta_2

alpha1 <- c(-Inf, 0, Inf)
alpha2 <- c(-Inf, -1, 0, 1, Inf)

Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
yhelp <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)

y1star <- mu1 + yhelp[,1]
y2star <- mu2 + yhelp[,2]

y1 <- cut(y1star, alpha1, labels = FALSE) - 1
y2 <- cut(y2star, alpha2, labels = FALSE) - 1

d <- data.frame(
  "x1" = x1,
  "x2" = x2,
  "x3" = x3,
  "x4" = x4,
  "x5" = x5,
  "x6" = x6,
  "y1" = y1,
  "y2" = y2
)

truth.coefs <- c(beta_1, rep(0,8), beta_2, rep(0,8), rho, 0, -1, 0, 1)

## Estimation in gamlss2.
## Response.
d$Y <- cbind(y1, y2)

f <- Y ~ la(x1,type=2) + la(x2,type=2) + la(x3,type=2) |
  la(x4,type=2) + la(x5,type=2) + la(x6,type=2)

source("fitting_routine.R")

## Estimate model.
if(!file.exists("simmodel.rds")) {
  b <- gamlss2(f, family = dk4, data = d, eps = .00001)
  saveRDS(b, file = "simmodel.rds")
} else {
  b <- readRDS("simmodel.rds")
}

pdf(file = "simulation_paths.pdf", width = 14 - 2*14/5, height = 9)

ci <- 0.7
lwd <- 2

par(mfrow = c(2, 3), mar = c(4, 4, 4, 3), oma = c(0, 4, 0, 0))

plot_lasso(b, which = "coefficients",
  term = "mu1.la(x1", spar = FALSE,
  zoom = c(8, 3), main = "x1",
  info = TRUE, cex.info = ci, lwd = lwd)

mtext(bquote(eta[1]), side = 2, line = 5, font = 2, cex = 1.5)

plot_lasso(b, which = "coefficients",
  term = "mu1.la(x2", spar = FALSE,
  zoom = c(3, 3), main = "x2",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu1.la(x3", spar = FALSE,
  zoom = c(3, 3), main = "x3",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(x4", spar = FALSE,
  zoom = c(8, 4), main = "x4",
  info = TRUE, cex.info = ci, lwd = lwd)

mtext(bquote(eta[2]), side = 2, line = 5, font = 2, cex = 1.5)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(x5", spar = FALSE,
  zoom = c(5, 3), main = "x5",
  info = TRUE, cex.info = ci, lwd = lwd)

plot_lasso(b, which = "coefficients",
  term = "mu2.la(x6", spar = FALSE,
  zoom = c(5, 3), main = "x6",
  info = TRUE, cex.info = ci, lwd = lwd)

dev.off()

