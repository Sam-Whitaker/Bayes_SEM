################################################################################
#                      Run inference code for Ebola data                       #
################################################################################

source("ebolaInf.r")

data <- read.csv("ebolaData.csv")$x

set.seed(1)
system.time(
  ebola_output <- ebola_storvik(
    parts = 4e6, dt = 0.1, a = c(2, 5, 10, 0.85, 5),
    b = c(50000, 4.6, 10, 0.75, 0.2)
  )
)

################################################################################
#                Run the Poisson leap code to produce forecasts                #
################################################################################

source("ebolaForecast.r")

set.seed(1)
ebola_forecast <- matrix(0, ncol = 5, nrow = 4e6)
for (i in 1:5) {
  print(i)
  paras <- ebola_output[[1]][, , i]
  state <- ebola_output[[2]][, , i]
  run <- t(
    simplify2array(
      mclapply(
        1:4e6, ebola_pol, state, paras[, 1:3], paras[, 4:5], 1, 0.1,
        mc.cores = detectCores() - 8
      )
    )
  )
  ebola_forecast[, i] <- rnbinom(
    4e6,
    mu = paras[, 4] * run,
    size = paras[, 5]
  )
}
