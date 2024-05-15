################################################################################
#                Run the Poisson leap code to produce forecasts                #
################################################################################

library(parallel)
library(boot)

source("forecast.R")

set.seed(1)
covid_forecast <- matrix(0, ncol = 5, nrow = 4e6)
for (i in 1:5){
  print(i)
  paras <- 
  state <-
  run <- t(
    simplify2array(
      mclapply(
        1:4e6, covid_pol, state, paras[, 1:4], paras [5:6], 1, 0.1, sepi,
        mc.cores = detectCores() - 8
      )
    )
  )
  covid_forecast[, i] <- rnbinom(
    4e6,
    mu = run[, 2] * run[, 1],
    size = 1 / paras[, 6]^2
  )
}