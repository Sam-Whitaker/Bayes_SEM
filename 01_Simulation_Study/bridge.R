################################################################################
#    Parallelised implementation of the bridge for one observation interval    #
################################################################################

ss_bridge = function(j, matx, mattheta, lbeta, obs, inter, dt, stoic) {
  # Input:
  #   - j: Indeced for particles (used for parallelisation) e.g. 1:100 for N=100
  #   - matx: [Nx2] matrix containing x_t for each particle
  #   - mattheta: [Nx3] matrix containing (gamma,lambda,rho) for each particle
  #   - lbeta: Vector (length N) containing log(beta_t) for each particle
  #   - obs: the next observation
  #   - inter: the inter-observation time
  #   - dt: the Poisson leap time step
  #   - stoic: the stoichiometry matrix
  len <- inter / dt + 1
  # Initialise by extracting correct components from each of the inputs
  x <- matx[j, ]
  theta <- mattheta[j, ]
  xmat <- matrix(0, ncol = 2, nrow = len)
  xmat[1, ] <- x
  # Initialise storage of sufficient statistics and log(beta)
  cumr <- matrix(0, ncol = 2, nrow = len)
  cumg <- matrix(0, ncol = 2, nrow = len)
  cumg[1, 1] <- xmat[1, 1] * xmat[1, 2]
  cumg[1, 2] <- xmat[1, 2]
  llr <- 0
  log_inf <- vector("numeric", len)
  log_inf[1] <- lbeta[j]
  sumsq <- 0
  # Begin loop
  for (i in 2:len) {
    # Calculate hazard and conditioned hazard
    h <- c(
      xmat[i - 1, 1] * xmat[i - 1, 2] * exp(log_inf[i - 1]),
      xmat[i - 1, 2] * theta[1]
    )
    hs <- cond_haz(
      xmat[i - 1, ], theta, log_inf[i - 1], obs, (i - 2) * dt, inter,
      cumr[i - 1, ]
    )
    # Simulate number of reactions
    rnew <- rpois(2, hs * dt)
    # Update state and sufficient statistics
    xmat[i, ] <- xmat[i - 1, ] + stoic %*% rnew
    cumr[i, ] <- cumr[i - 1, ] + rnew
    cumg[i, ] <- cumg[i - 1, ] + c(xmat[i, 1] * xmat[i, 2], xmat[i, 2])
    # Update log-likelihood ratio
    llr <- llr + sum(dpois(rnew, h * dt, log = TRUE)) -
      sum(dpois(rnew, hs * dt, log = TRUE))
    # Increment Brownian motion and calc sumsq
    log_inf[i] <- log_inf[i - 1] + (1 / sqrt(theta[2])) * rnorm(1, 0, sqrt(dt))
    sumsq <- sumsq + (log_inf[i - 1] - log_inf[i])^2
  }
  # Calculate the logged weight
  lwt <- dbinom(obs, cumr[len, 1], theta[3], log = TRUE) + llr
  return(c(xmat[len, ], cumr[len, ], cumg[len - 1, ], sumsq, log_inf[len], lwt))
}