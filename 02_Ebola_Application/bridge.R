################################################################################
#    Parallelised implementation of the bridge for one observation interval    #
################################################################################

ebola_bridge = function(i, state, theta, phi, obs, inter, dtau, sepi = stoic) {
  # Input:
  #  - i: the particle number index
  #  - state, theta, phi are all matrices with rows as in vectors above
  #  - obs: the next observation
  #  - inter: the inter-observation time
  #  - dtau: the Poisson leap time step
  #  - sepi: the stoichiometry matrix
  # Initialise storage
  len <- inter / dtau + 1
  btheta <- theta[i, ]
  bphi <- phi[i, ]
  xmat <- matrix(0, ncol = 3, nrow = len)
  xmat[1, ] <- state[i, ]
  sumr <- matrix(0, ncol = 3, nrow = len)
  sumg <- matrix(0, ncol = 3, nrow = len)
  sumg[1, ] <- c(xmat[1, 1] * xmat[1, 3], xmat[1, 2], xmat[1, 3])
  llr <- 0
  # Begin loop
  for (k in 1:(len - 1)) {
    # Calculate the hazard and conditioned hazard
    h <- hazard(xmat[k, ], btheta)
    hs <- cond_hazard(
      xmat[k, ], btheta, bphi, obs, (k - 1) * dtau, inter, sumr[k, ]
    )
    # Simulate incidence
    rnew <- rpois(3, hs * dtau)
    # Update state and summaries
    xmat[k + 1, ] <- xmat[k, ] + sepi %*% rnew
    if (any(xmat[k + 1, ] < 0) == TRUE) {
      xmat[k + 1, which(xmat[k + 1, ] < 0)] <- 0
      llr <- -Inf
    }
    sumr[k + 1, ] <- sumr[k, ] + rnew
    sumg[k + 1, ] <- sumg[k, ] +
      c(
        xmat[k + 1, 1] * xmat[k + 1, 3], xmat[k + 1, 2], xmat[k + 1, 3]
      )
    # Calculate the log-likelihood ratio
    llr <- llr + sum(dpois(rnew, h * dtau, log = TRUE)) -
      sum(dpois(rnew, hs * dtau, log = TRUE))
  }
  # Compute the log-weight
  lwt <- llr + dnbinom(
    obs, mu = phi[1] * sumr[len, 2], size = bphi[2], log = TRUE
  )
  return(c(xmat[len, ], sumr[len, ], sumg[(len - 1), ], lwt))
}