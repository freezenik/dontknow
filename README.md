# R package dontknow

`dontknow` provides a `gamlss2` family for models with an explicit "don't know"
response option, together with a simulation helper and an example data set.

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
f <- Y ~ -1 + s(x1) | -1 + s(x2) | 1
b <- gamlss2(f, family = DK(4), data = d)

summary(b)
```
