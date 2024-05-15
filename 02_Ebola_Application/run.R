################################################################################
#                      Run inference code for Ebola data                       #
################################################################################

library(parallel)
library(boot)

source("hazard.R")
source("bridge.R")
source("particle_filter.R")

data <- read.csv("data.csv")$x

set.seed(1)
system.time(
  ebola_output <- ebola_storvik(
    parts = 4e6, dt = 0.1, a = c(2, 5, 10, 0.85, 5),
    b = c(50000, 4.6, 10, 0.75, 0.2)
  )
)