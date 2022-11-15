# devtools::install_github("sooyongl/IMPS2022")
library(IMPS2022)
library(rstan)

# generate data ------------------------------------------
sdat <-
  makeDat(N = 500, R2Y = 0.2, R2eta = 0.5, omega = 0.5, tau0 = 0.3, tau1 = -0.15,
          lambda = 0.5, nsec = 100, lvmodel = '2pl')

stan_code <- "

"

fit <- rstan::stan(model_code = stan_code, data   = sdat$stan_dt)


