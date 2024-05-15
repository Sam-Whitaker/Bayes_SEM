################################################################################
#                   Run inference code for Simulation Study                    #
################################################################################

library(parallel)

source("hazard.R")
source("bridge.R")
source("particle_filter.R")

data <- read.csv("data.csv")$x

set.seed(1)
system.time(
  ss_output <- ss_storvik(
    parts = 1e4, dt = 0.1
  )
)