################################################################################

#                Functions for inference for Ebola application                 #

################################################################################

library(parallel)
library(boot)

################################################################################
#  Helper functions to generate draws from a multivariate Normal distribution  #
################################################################################

sqrtmat = function(v) {
  spec <- svd(v)
  return(spec$u %*% diag(sqrt(spec$d)) %*% t(spec$u))
}

rmvn = function(m, v) {
  p <- length(m)
  z <- rnorm(p)
  return(m + sqrtmat(v) %*% z)
}


################################################################################
#     Define hazard function, stoichiometry matrix, P and condition hazard     #
################################################################################

hazard = function(state, theta) {
  # state is a length-3 vector (S,E,I)
  # rate is a length-3 vector (beta,kappa,gamma)
  c(
    theta[1] * state[1] * state[3],
    theta[2] * state[2],
    theta[3] * state[3]
  )
}

stoic <- matrix(
  c(
    -1, 0, 0,
    1, -1, 0,
    0, 1, -1
  ), ncol = 3, nrow = 3, byrow = TRUE
)

cond_hazard = function(state, theta, phi, obs, curr, inter, rcurr) {
  # phi is a length-2 vector (rho,nu)
  h <- hazard(state, theta)
  mu <- phi[1] * (rcurr[2] + h[2] * (inter - curr))
  out <- c(
    h[1],
    h[2] + phi[1] * h[2] *
      (phi[1]^2 * h[2] * (inter - curr) + (mu + mu^2 / phi[2]))^(-1) *
      (obs - mu),
    h[3]
  )
  if (is.na(out[2])) {
    out[2] <- 0
  }
  out[out < 0] <- 0
  return(as.numeric(out))
}

################################################################################
#    Parallelised implementation of the bridge for one observation interval    #
################################################################################

ebola_bridge = function(i, state, theta, phi, obs, inter, dtau, sepi = stoic) {
  # state, theta, phi are all be matrices with rows as in vectors above
  len <- inter / dtau + 1
  btheta <- theta[i, ]
  bphi <- phi[i, ]
  xmat <- matrix(0, ncol = 3, nrow = len)
  xmat[1, ] <- state[i, ]
  sumr <- matrix(0, ncol = 3, nrow = len)
  sumg <- matrix(0, ncol = 3, nrow = len)
  sumg[1, ] <- c(xmat[1, 1] * xmat[1, 3], xmat[1, 2], xmat[1, 3])
  llr <- 0
  count <- 1
  for (k in 1:(len - 1)) {
    h <- hazard(xmat[k, ], btheta)
    hs <- cond_hazard(
      xmat[k, ], btheta, bphi, obs, (k - 1) * dtau, inter, sumr[k, ]
    )
    rnew <- rpois(3, hs * dtau)
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
    llr <- llr + sum(dpois(rnew, h * dtau, log = TRUE)) -
      sum(dpois(rnew, hs * dtau, log = TRUE))
    count <- count + 1
  }
  lwt <- llr + dnbinom(
    obs, mu = phi[1] * sumr[len, 2], size = bphi[2], log = TRUE
  )
  return(c(xmat[len, ], sumr[len, ], sumg[(len - 1), ], lwt))
}

################################################################################
#                  Particle filter code for Ebola application                  #
################################################################################

ebola_storvik = function(
  parts = 100000, dat = data, x0 = c(44326, 15, 10), inter = 1, dt = 1 / 7,
  a = c(0.2, 5, 10, 0.85, 5), b = c(5000, 4.6, 10, 0.75, 0.1), sepi = stoic
) {
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
