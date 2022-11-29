
# IMPS 2022 Proceedings

## Overview

This package is for generating simulation data for FLPS model which is
introduced in **IMPS 2022 proceedings**. More information can be found
in **Fully Latent Principal Stratification: combining PS with
model-based measurement models**.

## Install

``` r
devtools::install_github("sooyongl/IMPS2022")
```

## Generate Data

``` r
library(IMPS2022)
library(rstan)

sdat <- makeDat(
  N = 200,        # number of total sample (Half is for the treatment group)
  R2Y = 0.2,      # Proportion of outcome explained by covariates
  R2eta = 0.5,    # Proportion of latent variance explained by covariates
  omega = 0.5,    # Effect of latent variable on outcome
  tau0 = 0.3,     # Difference in outcome between Treatment and Control
  tau1 = -0.15,   # Principal effects
  lambda = 0.5,   # Missingness on item responses
  nsec = 20,      # Number of items
  lvmodel = '2pl' # Latent variable models
  )
```

## Run Stan for FLPS

``` r
# You can load the raw stan script corresponding to the type of measurement models from the package
stan_model <- mk_stanmodel("2pl") # "2pl", "grm", "gpcm" available

# Run stan based on the generated data and the FLPS stan script
fit <- rstan::stan(model_code = stan_code, data = sdat$stan_dt)

summary(fit)
```
