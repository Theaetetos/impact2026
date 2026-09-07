# impact2026

Worked examples of (controlled) interrupted time series, difference-in-differences,
and propensity score matching with G-computation for valuing a managed care
intervention on a synthetic dataset. Supporting material for the SOA ImpACT 2026
session.

## Setup

```r
install.packages("renv")
renv::init()          # creates renv.lock, renv/, and .Rprofile
devtools::load_all()
devtools::test()
```

## Layout

```
R/                 package functions
man/               roxygen-generated documentation
tests/testthat/    testthat suite (edition 3)
DESCRIPTION        metadata + dependencies (renv reads Imports/Suggests from here)
NAMESPACE          roxygen-generated exports
```
