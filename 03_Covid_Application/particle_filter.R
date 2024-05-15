################################################################################
#                  Particle filter code for Covid application                  #
################################################################################

covid_storvik = function(
    parts = 100000, dat = data, x0, inter = 1, dt = 1 / 7, 
    a = c(11.088, 4, 4, 3, 0), b = c(2.192, 1, 1, 2, 0.5), sepi = stoic
    ){
  # Input:
  #  - parts: the number of particles to use
  #  - dat: the data
  #  - x0: the initial state (as a matrix)
  #  - inter: the inter-observation time
  #  - dt: the Poisson leap time step
  #  - a, b: prior hyper-parameters
  #  - sepi: the stoichiometry matrix
  # Initialise with prior draws
  theta <- matrix(0, ncol = 7, nrow = parts)
  # Theta contains lbeta, gamma, lambda_beta, lambda_rho, rho_t, nu, and R0
  theta[, 1] <- rep(0, N)
  theta[, 2] <- rgamma(N, shape = a[1], rate = b[1])
  theta[, 3] <- rgamma(N, shape = a[2], rate = b[2])
  theta[, 4] <- rgamma(N, shape = a[3], rate = b[3])
  theta[, 5] <- rbeta(N, a[4], b[4])
  theta[, 6] <- runif(N, a[5], b[5])
  theta[, 7] <- exp(theta[, 1]) / theta[, 2]
  xmat <- x0
  # Set up storage of parameters and states for forecasting
  cstore <- array(dim = c(parts, 7, 5))
  xstore <- array(dim = c(parts, 2, 5))
  # Set up storage for summaries
  # 9 means followed by 9 lower quantiles followed by 9 upper quantiles
  sum_traj <- matrix(0, ncol = 27, nrow = (length(dat) + 1))
  sum_traj[1, ] <- c(
    apply(cbind(theta, xmat), 2, mean),
    c(t(apply(cbind(theta, xmat), 2, quantile, c(0.025, 0.975))))
  )
  lat_traj <- matrix(0, ncol = length(data), nrow = parts)
  obs_traj <- matrix(0, ncol = length(data), nrow = parts)
  sumr <- matrix(0, ncol = 2, nrow = N)
  sumg <- matrix(0, ncol = 2, nrow = N)
  bsumsq <- rep(0,N)
  rsumsq <- rep(0,N)
  #wts = vector("numeric",N)
  # Calculate Liu-West h and a
  hjit <- sqrt(1 - ((3 * 0.99 - 1) / (2 * 0.99))^2)
  ajit <- sqrt(1 - hjit^2)
  # Loop over all observations
  for (i in 1:length(data)) {
    message(i)
    # Liu-West update for nu
    lnu <- log(theta[, 6])
    m <- ajit * lnu + (1 - ajit) * mean(lnu)
    V <- sum((lnu - mean(lnu))^2) / N
    theta[, 6] <- exp(rnorm(N, m, hjit * sqrt(V)))
    # Propagate dynamic processes
    prop <- t(
      simplify2array(
        mclapply(
          1:parts, covid_bridge, state = xmat, theta = theta[, 1:4], 
          phi = theta[, 5:6], obs = dat[i], inter = inter, dtau = dt, 
          sepi = sepi, mc.cores = detectCores() - 8
        )
      )
    )
    # Sort output
    xmat <- prop[, 1:2]
    sumr <- sumr + prop[, 3:4]
    sumg <- sumg + prop[, 5:6]
    bsumsq <- bsumsq + prop[, 7]
    theta[, 1] <- prop[, 8]
    wts <- exp(prop[, 9])
    theta[, 5] <- prop[, 10]
    rsumsq <- rsumsq + prop[, 11]
    # Resample and update summaries
    wts[is.na(wts)] <- 0
    pick <- sample(1:parts, parts, replace = TRUE, prob = wts)
    xmat <- xmat[pick, ]
    sumr <- sumr[pick, ]
    sumg <- sumg[pick, ]
    bsumsq <- bsumsq[pick]
    rsumsq <- rsumsq[pick]
    theta <- theta[pick, ]
    # Rejuvenate theta
    theta[, 2] <- rgamma(
      parts, shape = a[1] + sumr[, 2], rate = b[1] + dtau * sumg[, 2]
    )
    theta[, 3] <- rgamma(
      parts, shape = a[2] + 0.5 * i * inter / dtau, 
      rate = b[2] + 0.5 / dtau * bsumsq
    )
    theta[, 4] <- rgamma(
      parts, shape = a[3] + 0.5 * i * inter / dtau, 
      rate = b[3] + 0.5 / dtau * rsumsq
    )
    theta[, 7] <- exp(theta[, 1]) / theta[, 2]
    sum_traj[(i + 1), ] <- c(
      apply(cbind(theta, xmat), 2, mean),
      apply(cbind(theta, xmat), 2, quantile, 0.025),
      apply(cbind(theta, xmat), 2, quantile, 0.975)
    )
    lat_traj[, i] <- prop[pick, 3]
    obs_traj[, i] <- rnbinom(
      parts, mu = theta[, 5] * lat_traj[, i], size = 1 / (theta[, 6]^2))
    # Store paras and states for last five observations for forecasting
    if (any(20:24 == i)) {
      cstore[, , i - 19] <- theta
      xstore[, , i - 19] <- xmat
    }
  }
  return(list(theta = theta, trajs = sum_traj, lat_traj, obs_traj, cstore, xstore))
}