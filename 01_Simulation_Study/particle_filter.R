################################################################################
#                     Particle filter for simulation study                     #
################################################################################

ss_storvik = function(
    parts = 10000, dat = data, x0 = c(762, 5), inter = 1, dt = 0.01, 
    a = c(-6.5, 11, 15, 90), b = c(0.5, 20, 0.14, 15), sepi = stoic
) {
  # Input:
  #   - parts: the number of particles to use
  #   - dat: the simulated data
  #   - x0: the initial state
  #   - inter: the inter-observation time
  #   - dt: the Poisson leap time step
  #  - a, b: prior hyper-parameters
  #   - sepi: the stoichiometry matrix
  # Initialise with prior draws
  theta <- matrix(0, nrow = parts, ncol = 4)
  # theta contains beta, gamma, lambda, rho
  theta[, 1] <- rnorm(N, a[1], b[1])
  theta[, 2] <- rgamma(parts, shape = a[2], rate = b[2])
  theta[, 3] <- rgamma(parts, shape = a[3], rate = b[3])
  theta[, 4] <- rbeta(parts, a[4], b[4])
  xmat <- matrix(x0, nrow = parts, ncol = 2, byrow = TRUE)
  # Set up storage for summaries
  # 6 means followed by 6 lower quantiles followed by 6 upper quantiles
  sum_traj <- matrix(0, ncol = 18, nrow = (length(data) + 1))
  sum_traj[1, ] <- c(
    apply(cbind(theta, xmat), 2, mean),
    c(t(apply(cbind(theta, xmat), 2, quantile, c(0.025, 0.975))))
  )
  lat_traj <- matrix(0, ncol = 10, nrow = parts)
  obs_traj <- matrix(0, ncol = 10, nrow = parts)
  sumr <- matrix(0, ncol = 2, nrow = parts)
  sumg <- matrix(0, ncol = 2, nrow = parts)
  sumsq <- rep(0, parts)
  # Loop over all observations
  for (i in seq_along(dat)) {
    message(i)
    # Propagate dynamic processes
    prop <- t(
      simplify2array(
        mclapply(
          1:parts, ss_bridge, xmat, theta[, 2:4], theta[, 1], dat[i], inter,
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
      parts, shape = a[2] + sumr[, 2], rate = b[2] + dt * sumg[, 2]
    )
    theta[, 3] <- rgamma(
      parts, shape = a[3] + 0.5 * i * inter / dt, rate = b[3] + 0.5 / dt * sumsq
    )
    theta[, 4] <- rbeta(
      parts, a[4] + sum(data[1:i]), b[4] + sumr[, 1] - sum(data[1:i])
    )
    sum_traj[i + 1, ] <- c(
      apply(cbind(theta, xmat), 2, mean),
      apply(cbind(theta, xmat), 2, quantile, 0.025),
      apply(cbind(theta, xmat), 2, quantile, 0.975)
    )
    lat_traj[, i] <- prop[pick, 3]
    obs_traj[, i] <- rbinom(parts, lat_traj[, i], theta[, 4])
  }
  return(list(theta = theta[, 2:4], trajs = sum_traj, lat_traj, obs_traj))
}
