################################################################################
#    Parallelised implementation of the bridge for one observation interval    #
################################################################################

covid_bridge = function(j, state, theta, phi, obs, inter, dtau, sepi = stoic) {
  # Input:
  #  - j: the particle number index
  #  - state, theta, phi are all matrices with rows as in vectors above
  #  - obs: the next observation
  #  - inter: the inter-observation time
  #  - dtau: the Poisson leap time step
  #  - sepi: the stoichiometry matrix
  # Initialise storage
  len <- inter / dtau + 1
  btheta <- theta[j, ]
  bphi <- phi[j, ]
  xmat <-  matrix(0, ncol = 2, nrow = len)
  xmat[1, ] <- state[j, ]
  sumr <- matrix(0, ncol = 2, nrow = len)
  sumg <- matrix(0, ncol = 2, nrow = len)
  sumg[1, ] <- c(xmat[1, 1] * xmat[1, 2], xmat[1, 2])
  llr <- 0
  bsumsq <- 0
  rsumsq <- 0
  # Begin loop
  for (i in 1:(len - 1)){
    # Store the current beta and rho for later calculation
    linfold <- btheta[1]
    lrold <- logit(bphi[1])
    # Calculate the hazard and conditioned hazard
    h <- c(
      xmat[i, 1] * xmat[i, 2] * exp(btheta[1]) / 19453561,
      xmat[i, 2] * theta[2]
    )
    hs <- cond_hazard(
      xmat[i, ], btheta, bphi, obs, (i - 1) * dtau, inter, sumr[i, ]
    )
    # Simulate incidence
    rnew <- rpois(2, hs * dtau)
    # Update state and summaries
    xmat[i + 1, ] <- xmat[i, ] + sepi %*% rnew
    if (any(xmat[i + 1, ] < 0) == TRUE) {
      xmat[i + 1, which(xmat[i + 1, ] < 0)] <- 0
      llr <- -Inf
    }
    sumr[i + 1, ] <- sumr[i, ] + rnew
    sumg[i + 1, ] <- sumg[i, ] + 
      c(
        xmat[i + 1, 1] * xmat[i + 1, 2], xmat[i + 1, 2]
      )
    
    # Calculate the log-likelihood ratio
    llr <- llr + sum(dpois(rnew, h * dtau, log = TRUE)) -
      sum(dpois(rnew, hs * dtau, log = TRUE))
    # Increment Brownian motion and calculate the sum of squares
    btheta[1] <- rnorm(1, btheta[1], sqrt(dtau / btheta[3]))
    bphi[1] <- inv.logit(rnorm(1, logit(bphi[1]), sqrt(dtau / btheta[4])))
    bsumsq <- bsumsq + (linfold - btheta[1])^2
    rsumsq <- rsumsq + (lrold - logit(bphi[1]))^2
  }
  # Compute the log-weight
  lwt <- dnbinom(
    obs, mu = bphi[1] * sumr[len, 1], size = 1 / (bphi[2]^2), log = TRUE
  ) + llr
  return(
    c(
      xmat[len, ], sumr[len, ], sumg[len - 1, ], bsumsq, btheta[1], lwt,
      bphi[1], rsumsq
    )
  )
}
