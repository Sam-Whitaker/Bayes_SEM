################################################################################
#                      Run inference code for Covid data                       #
################################################################################

library(parallel)
library(boot)

source("hazard.R")
source("bridge.R")
source("particle_filter.R")

data <- read.csv("data.csv")$x

set.seed(1)
parts = 4e6
I0 = runif(parts,-16,-9)
R0 = rnorm(parts,0.5,0.15)
while (any(R0<0 | R0>1) == TRUE){
  R0[R0<0|R0>1] = rnorm(
    length(R0[R0<0|R0>1]),0.5,0.15
  )
}
N = 19453561
R = R0*N
I = exp(I0+log(N-R))
S = N - R - I
x0 = cbind(S,I)

system.time(
  covid_output <- covid_storvik(
    parts = 4e6, dt = 0.1, a = c(2, 5, 10, 0.85, 5),
    b = c(50000, 4.6, 10, 0.75, 0.2)
  )
)