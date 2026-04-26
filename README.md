# R package dontknow

`dontknow` provides a `gamlss2` family for regression models with an explicit
"don't know" response option. The response is represented as
`Y = cbind(y1, y2)`, where `y1` indicates whether the DK option was chosen and
`y2` contains the ordinal response among non-DK observations.

`DK(k)` fits a joint latent bivariate normal model for:
- the DK decision
- the ordinal outcome
- their dependence through a correlation parameter

The package also contains `sim_DK()` for generating example data from the same
model structure.

The package depends on `gamlss2`, which can be installed from R-universe:

```r
install.packages(
  "gamlss2",
  repos = c(
    "https://gamlss-dev.R-universe.dev",
    "https://cloud.R-project.org"
  )
)
```

The package itself can then be installed from GitHub:

```r
remotes::install_github("freezenik/dontknow")
```

Minimal example:

```r
library(gamlss2)
library(dontknow)

d <- sim_DK(500)

# y1: DK indicator, y2: ordinal category among non-DK responses
head(d$Y)

# remove intercepts in the latent mean predictors
f <- Y ~ -1 + s(x1) | -1 + s(x2) | 1
b <- gamlss2(f, family = DK(4), data = d)

summary(b)

# fitted distribution parameters
par <- predict(b, type = "parameter")
head(par$mu1)
head(par$mu2)
head(par$rho)

# ordered thresholds on the latent scale
head(family(b)$alpha(par))
```

Implementation note:

`DK(4)` uses compiled code by default for the log-likelihood, score, and
observed Hessian. When score and Hessian are both enabled, the family also uses
a fused compiled `z_weights` path for faster RS updates in `gamlss2`.
