################################################################################
#                  Particle filter code for Ebola application                  #
################################################################################

ebola_storvik = function(
    parts = 100000, dat = data, x0 = c(44326, 15, 10), inter = 1, dt = 1 / 7,
    a = c(0.2, 5, 10, 0.85, 5), b = c(5000, 4.6, 10, 0.75, 0.1), sepi = stoic
) {
  # Input:
  #  - parts: the number of particles to use
  #  - dat: the data
  #  - x0: the initial state
  #  - inter: the inter-observation time
  #  - dt: the Poisson leap time step
  #  - a, b: prior hyper-parameters
  #  - sepi: the stoichiometry matrix
  #Initialise with prior draws
  cmat <- matrix(0, ncol = 5, nrow = parts)
  # cmat contains beta, omega, gamma, rho, nu
  cmat[, 1] <- rgamma(parts, shape = a[1], rate = b[1])
  cmat[, 2] <- rgamma(parts, shape = a[2], rate = b[2])
  cmat[, 3] <- rgamma(parts, shape = a[3], rate = b[3])
  cmat[, 4] <- rnorm(parts, a[4], b[4])
  cmat[, 4] <- inv.logit(cmat[, 4])
  cmat[, 5] <- rgamma(parts, shape = 5, rate = 0.2)
  xmat <- matrix(x0, ncol = 3, nrow = parts, byrow = TRUE)
  # Set up storage of parameters and states for forecasting
  cstore <- array(dim = c(parts, 5, 5))
  xstore <- array(dim = c(parts, 3, 5))
  # Set up storage for summaries
  # 8 means followed by 8 lower quantiles followed by 8 upper quantiles
  sum_traj <- matrix(0, ncol = 24, nrow = (length(dat) + 1))
  sum_traj[1, ] <- c(
    apply(cbind(cmat, xmat), 2, mean),
    c(t(apply(cbind(cmat, xmat), 2, quantile, c(0.025, 0.975))))
  )
  lat_traj <- matrix(0, ncol = length(dat), nrow = parts)
  obs_traj <- matrix(0, ncol = length(dat), nrow = parts)
  sumr <- matrix(0, ncol = 3, nrow = parts)
  sumg <- matrix(0, ncol = 3, nrow = parts)
  # Calculate Liu-West h and a
  hjit <- sqrt(1 - ((3 * 0.99 - 1) / (2 * 0.99))^2)
  ajit <- sqrt(1 - hjit^2)
  # Loop over all observations
  for (i in seq_along(dat)) {
    message(i)
    # Liu-West update for rho and nu
    lrho <- logit(cmat[, 4])
    lnu <- log(cmat[, 5])
    phi <- rbind(lrho, lnu)
    phibar <- apply(phi, 1, mean)
    phitild <- phi - phibar
    v <- 1 / parts * phitild %*% t(phitild)
    for (k in 1:parts) {
      m <- ajit * phi[, k] + (1 - ajit) * phibar
      draw <- rmvn(m, hjit^2 * v)
      lrho[k] <- draw[1]
      lnu[k] <- draw[2]
    }
    cmat[, 4] <- inv.logit(lrho)
    cmat[, 5] <- exp(lnu)
    # Propagate dynamic processes
    prop <- t(
      simplify2array(
        mclapply(
          1:parts, ebola_bridge, state = xmat, theta = cmat[, 1:3],
          phi = cmat[, 4:5], obs = dat[i], inter = inter, dtau = dt,
          sepi = sepi, mc.cores = detectCores() - 8
        )
      )
    )
    # Sort output
    xmat <- prop[, 1:3]
    sumr <- sumr + prop[, 4:6]
    sumg <- sumg + prop[, 7:9]
    wts <- exp(prop[, 10])
    # Resample and update summaries
    wts[is.na(wts)] <- 0
    pick <- sample(1:parts, parts, replace = TRUE, prob = wts)
    xmat <- xmat[pick, ]
    sumr <- sumr[pick, ]
    sumg <- sumg[pick, ]
    cmat <- cmat[pick, ]
    # Rejuvenate theta
    cmat[, 1] <- rgamma(
      parts, shape = a[1] + sumr[, 1], rate = b[1] + dt * sumg[, 1]
    )
    cmat[, 2] <- rgamma(
      parts, shape = a[2] + sumr[, 2], rate = b[2] + dt * sumg[, 2]
    )
    cmat[, 3] <- rgamma(
      parts, shape = a[3] + sumr[, 3], rate = b[3] + dt * sumg[, 3]
    )
    sum_traj[(i + 1), ] <- c(
      apply(cbind(cmat, xmat), 2, mean),
      apply(cbind(cmat, xmat), 2, quantile, 0.025),
      apply(cbind(cmat, xmat), 2, quantile, 0.975)
    )
    lat_traj[, i] <- prop[pick, 5]
    obs_traj[, i] <- rnbinom(
      parts, mu = cmat[, 4] * lat_traj[, i], size = cmat[, 5]
    )
    # Store paras and states if at time 45, 46, 47, 48 or 49 (for forecasting)
    if (i %in% 45:49) {
      cstore[, , i - 44] <- cmat
      xstore[, , i - 44] <- xmat
    }
  }
  return(list(theta = cmat, trajs = sum_traj, lat_traj, obs_traj))
}