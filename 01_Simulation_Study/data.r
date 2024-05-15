################################################################################

#             Code to produce synthetic data for simulation study              #

################################################################################


################################################################################
#        The Poisson leap for SIR model with time-varying infection rate       #
################################################################################

poisson_leap = function(x0, theta, lbeta0, max_time, delta_t) {
  # Input:
  #   - x0: the initial state of the system
  #   - theta: the model parameters (gamma,lambda)
  #   - lbeta0: the initial value of the infection rate process
  #   - max_time: the maximum time to run the algorithm
  #   - delta_t: the time step to be used
  len <- max_time / delta_t + 1
  # Initialise storage for x, n, t and lbeta
  x <- matrix(0, nrow = len, ncol = 2)
  x[1, ] <- x0
  rmat <- matrix(0, nrow = len, ncol = 2)
  time <- vector("numeric")
  time[1] <- 0
  lbeta <- vector("numeric")
  lbeta[1] <- lbeta0
  # Define the Stoichiometry matrix
  stoic <- matrix(
    c(-1, 0, 1, -1),
    ncol = 2, nrow = 2, byrow = TRUE
  )
  # Start the counter
  count <- 1
  # Begin loop
  while (time[count] < max_time) {
    # Calculate the hazard function
    h <- c(
      exp(lbeta[count]) * x[count, 1] * x[count, 2],
      theta[1] * x[count, 2]
    )
    h0 <- sum(h)
    # Check for absorbing state
    if (h0 <= 0) {
      break
    }
    # Simulate number of new reactions over the interval
    rnew <- rpois(2, h * delta_t)
    # Update cumulative incidence and state
    rmat[count + 1, ] <- rmat[count, ] + rnew
    x[count + 1, ] <- x[count, ] + stoic %*% rnew
    # Check for negativity of state
    if (any(x[count + 1, ] < 0) == TRUE) {
      x[count + 1, ][x[count + 1, ] < 0] <- 0
    }
    time[count + 1] <- time[count] + delta_t
    lbeta[count + 1] <- rnorm(
      1, lbeta[count], 1 / sqrt(theta[2]) * sqrt(delta_t)
    )
    count <- count + 1
  }
  out <- data.frame(time, x[1:count, ], lbeta, rmat[1:count, ])
  names(out) <- c("time", "S", "I", "lbeta", "r1", "r2")
  return(out)
}

################################################################################
#       Run with x0=(762,5), theta=(0.5,100), log(beta0)=-6 and dt=0.001       #
################################################################################

set.seed(1)
dat <- poisson_leap(c(762, 5), c(0.5, 100), -6, 30, 0.001)

################################################################################
#                    Extract the log(beta) and n1 processes                    #
################################################################################

lbeta <- dat$lbeta[seq(1, length(dat$r1), 1000)]
(data <- dat$r1[seq(1, length(dat$r1), 1000)])

################################################################################
#          Subject latent processes to binomial error with rho = 0.9           #
################################################################################

y <- rbinom(10, data[2:11] - data[1:10], 0.9)
