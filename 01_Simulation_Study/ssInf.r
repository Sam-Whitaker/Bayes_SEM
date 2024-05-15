################################################################################

#                 Functions for inference for simulation study                 #

################################################################################

################################################################################
#                           The stoichiometry matrix                           #
################################################################################

stoic <- matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE)

################################################################################
#                        The conditioned hazard function                       #
################################################################################

cond_haz = function(x, theta, lbeta, obs, curr, inter, rcur) {
  pmat <- matrix(c(1, 0), ncol = 2, nrow = 1)
  h <- c(x[1] * x[2] * exp(lbeta), x[2] * theta[1])
  out <- h + theta[3] * diag(h) %*% t(pmat) %*%
    (
      theta[3]^2 %*% pmat %*% diag(h) %*% t(pmat) * (inter - curr) +
        theta[3] * (1 - theta[3]) * pmat %*% (rcur + h * (inter - curr))
    )^(-1) %*% (obs - theta[3] * pmat %*% (rcur + h * (inter - curr)))
  out[out <= 0] <- 1e-6
  out[out > 1e6] <- 1e6
  return(out)
}

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

################################################################################
#                    Particle filter for simulation example                    #
################################################################################

ss_storvik = function(
  data = y, parts = 10000, x0 = c(762, 5), inter = 1, dt = 0.01, sepi = stoic
) {
  # Input:
  #   - data: the simulated data
  #   - parts: the number of particles to use
  #   - x0: the initial state
  #   - inter: the inter-observation time
  #   - dt: the Poisson leap time step
  #   - sepi: the stoichiometry matrix
  
  # Initialise the storage
  theta <- matrix(0, nrow = parts, ncol = 4)
  xmat <- matrix(x0, nrow = parts, ncol = 2, byrow = TRUE)
  r1_traj <- matrix(0, ncol = 10, nrow = parts)
  obs_traj <- matrix(0, ncol = 10, nrow = parts)
  wts <- vector("numeric", parts)
  sumr <- matrix(0, ncol = 2, nrow = parts)
  sumg <- matrix(0, ncol = 2, nrow = parts)
  sumsq <- rep(0, parts)
  
  # Draw elements of theta from the prior
  theta[, 1] <- rep(-6, parts)
  theta[, 2] <- rgamma(parts, shape = 200, rate = 400)
  theta[, 3] <- rgamma(parts, shape = 1e4, rate = 1e2)
  theta[, 4] <- rbeta(parts, 81, 9)
  
  # Mean and quantile trajectories stored as matrix
  # In order: beta, gamma, rho, tau, S, I
  # 6 means followed by 6 lower quantiles followed by 6 upper quantiles
  sum_traj <- matrix(0, ncol = 18, nrow = (length(data) + 1))
  sum_traj[1, ] <- c(
    apply(cbind(theta, xmat), 2, mean),
    c(t(apply(cbind(theta, xmat), 2, quantile, c(0.025, 0.975))))
  )

  # Loop over all observations
  for (i in seq_along(data)) {
    message(i)
    # Propagate dynamic processes
    prop <- t(
      simplify2array(
        mclapply(
          1:parts, ss_bridge, xmat, theta[, 2:4], theta[, 1], data[i], inter,
          dt, sepi, mc.cores = detectCores() - 8
        )
      )
    )
    # Sort output
    xmat <- prop[, 1:2]
    sumr <- sumr + prop[, 3:4]
    sumg <- sumg + prop[, 5:6]
    sumsq <- sumsq + prop[, 7]
    theta[, 1] <- prop[, 8]
    wts <- exp(prop[, 9])
    # Resample and update summaries
    wts[is.na(wts)] <- 0
    pick <- sample(1:parts, parts, replace = TRUE, prob = wts)
    xmat <- xmat[pick, ]
    sumr <- sumr[pick, ]
    sumg <- sumg[pick, ]
    sumsq <- sumsq[pick]
    theta[, 1] <- theta[pick, 1]
    # Sample theta, lambda and phi
    theta[, 2] <- rgamma(
      parts, shape = 200 + sumr[, 2], rate = 400 + dt * sumg[, 2]
    )
    theta[, 3] <- rgamma(
      parts, shape = 1e4 + 0.5 * i * inter / dt, rate = 1e2 + 0.5 / dt * sumsq
    )
    theta[, 4] <- rbeta(
      parts, 81 + sum(data[1:i]), 9 + sumr[, 1] - sum(data[1:i])
    )
    sum_traj[i + 1, ] <- c(
      apply(cbind(theta, xmat), 2, mean),
      c(t(apply(cbind(theta, xmat), 2, quantile, c(0.025, 0.975))))
    )
    r1_traj[, i] <- prop[pick, 3]
    obs_traj[, i] <- rbinom(parts, r1_traj[, i], theta[, 4])
  }
  return(list(theta = theta[, 2:4], trajs = sum_traj, r1_traj, obs_traj))
}

################################################################################
# Run code for simulated data
################################################################################

library(parallel)
set.seed(1)
system.time(
  ss_output <- ssStorvik(N = 300000)
)
