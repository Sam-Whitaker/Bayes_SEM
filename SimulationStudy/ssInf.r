################################################################################

#                 Functions for inference for simulation study                 #

################################################################################

################################################################################
#                           The stoichiometry matrix                           #
################################################################################

Sepi <- matrix(c(-1, 0, 1, -1), nrow=2, byrow=TRUE)

################################################################################
#                        The conditioned hazard function                       #
################################################################################

cond_haz = function(x, theta, lbeta, yT, Tcur, TT, rcur) {
  Pt <- matrix(c(1, 0), ncol=2, nrow=1)
  h <- c(x[1]*x[2]*exp(lbeta), x[2]*theta[1])
  out <- h + theta[3]*diag(h)%*%t(Pt)%*%
    (
      theta[3]^2%*%Pt%*%diag(h)%*%t(Pt)*(TT-Tcur)+
        theta[3]*(1-theta[3])*Pt%*%(rcur+h*(TT-Tcur))
    )^(-1)%*%(yT-theta[3]*Pt%*%(rcur+h*(TT-Tcur)))
  out[out<=0] <- 1e-6
  out[out>1e6] <- 1e6
  return(out)
}

################################################################################
#    Parallelised implementation of the bridge for one observation interval    #
################################################################################

ssBridge = function(j, matx, mattheta, lbeta, yT, TT, dt, S) {
  # Input:
  #   - j: Indeced for particles (used for parallelisation) e.g. 1:100 for N=100
  #   - matx: [Nx2] matrix containing x_t for each particle
  #   - mattheta: [Nx3] matrix containing (gamma,lambda,rho) for each particle
  #   - lbeta: Vector (length N) containing log(beta_t) for each particle
  #   - yT: the next observation
  #   - TT: the inter-observation time
  #   - dt: the Poisson leap time step
  #   - S: the stoichiometry matrix
  N <- TT/dt + 1
  # Initialise by extracting correct components from each of the inputs
  x <- matx[j, ]
  theta <- mattheta[j, ]
  xmat <- matrix(0, ncol=2, nrow=N)
  xmat[1, ] <- x
  # Initialise storage of sufficient statistics and log(beta)
  cumr <- matrix(0, ncol=2, nrow=N)
  cumg <- matrix(0, ncol=2, nrow=N)
  cumg[1, 1] <- xmat[1, 1]*xmat[1, 2]
  cumg[1, 2] <- xmat[1, 2]
  llr <- 0
  log_inf <- vector("numeric", N)
  log_inf[1] <- lbeta[j]
  sumsq <- 0
  # Begin loop
  for (i in 2:N) {
    # Calculate hazard and conditioned hazard
    h <- c(xmat[i-1, 1]*xmat[i-1, 2]*exp(log_inf[i-1]), xmat[i-1, 2]*theta[1])
    hs <- cond_haz(
      xmat[i-1, ], theta, log_inf[i-1], yT, (i-2)*dt, TT, cumr[i-1, ]
    )
    # Simulate number of reactions
    rnew <- rpois(2, hs*dt)
    # Update state and sufficient statistics
    xmat[i, ] <- xmat[i-1, ] + S%*%rnew
    cumr[i, ] <- cumr[i-1, ] + rnew
    cumg[i, ] <- cumg[i-1, ] + c(xmat[i, 1]*xmat[i, 2], xmat[i, 2])
    # Update log-likelihood ratio
    llr <- llr + sum(dpois(rnew, h*dt, log=TRUE)) - 
      sum(dpois(rnew, hs*dt, log=TRUE))
    # Increment Brownian motion and calc sumsq
    log_inf[i] <- log_inf[i-1] + (1/sqrt(theta[2]))*rnorm(1, 0, sqrt(dt))
    sumsq <- sumsq + (log_inf[i-1] - log_inf[i])^2
  }
  # Calculate the logged weight
  lwt <- dbinom(yT, cumr[N, 1], theta[3], log=TRUE) + llr
  return(c(xmat[N, ], cumr[N, ], cumg[N-1, ], sumsq, log_inf[N], lwt))
}

################################################################################
#                    Particle filter for simulation example                    #
################################################################################

ssStorvik = function(data=y, N=10000, x0=c(762, 5), TT=1, dt=0.01, S=Sepi) {
  # Input:
  #   - data: the simulated data
  #   - N: the number of particles to use
  #   - x0: the initial state
  #   - TT: the inter-observation time
  #   - dt: the Poisson leap time step
  #   - S: the stoichiometry matrix
  
  # Initialise the storage
  theta <- matrix(0, nrow=N, ncol=4)
  xmat <- matrix(x0, nrow=N, ncol=2, byrow=TRUE)
  r1_traj <- matrix(0, ncol=10, nrow=N)
  obs_traj <- matrix(0, ncol=10, nrow=N)
  wts <- vector("numeric", N)
  sumr <- matrix(0, ncol=2, nrow=N)
  sumg <- matrix(0, ncol=2, nrow=N)
  sumsq <- rep(0, N)
  
  # Draw elements of theta from the prior
  theta[, 1] <- rep(-6, N)
  theta[, 2] <- rgamma(N, shape=200, rate=400)
  theta[, 3] <- rgamma(N, shape=1e4, rate=1e2)
  theta[, 4] <- rbeta(N, 81, 9)
  
  # Mean and quantile trajectories stored as matrix
  # In order: beta, gamma, rho, tau, S, I
  # 6 means followed by 6 lower quantiles followed by 6 upper quantiles
  sum_traj <- matrix(0, ncol=18, nrow=(length(data)+1))
  sum_traj[1, ] <- c(
    apply(cbind(theta, xmat), 2, mean),
    c(t(apply(cbind(theta, xmat), 2, quantile, c(0.025, 0.975))))
  )

  # Loop over all observations
  for (i in 1:length(data)) {
    message(i)
    # Propagate dynamic processes
    prop <- simplify2array(
      t(
        mclapply(
          1:N, ssBridge, xmat, theta[, 2:4], theta[, 1], data[i], TT, dt, S,
          mc.cores=detectCores()-8
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
    pick <- sample(1:N, N, replace=TRUE, prob=wts)
    xmat <- xmat[pick, ]
    sumr <- sumr[pick, ]
    sumg <- sumg[pick, ]
    sumsq <- sumsq[pick]
    theta[, 1] <- theta[pick, 1]
    # Sample theta, lambda and phi
    theta[, 2] <- rgamma(N, shape=200+sumr[, 2], rate=400+dt*sumg[, 2])
    theta[, 3] <- rgamma(N, shape=1e4+0.5*i*TT/dt, rate=1e2+0.5/dt*sumsq)
    theta[, 4] <- rbeta(N, 81+sum(data[1:i]), 9+sumr[, 1]-sum(data[1:i]))
    sum_traj[i+1, ] <- c(
      apply(cbind(theta, xmat), 2, mean),
      c(t(apply(cbind(theta, xmat), 2, quantile, c(0.025, 0.975))))
    )
    r1_traj[, i] <- prop[pick, 3]
    obs_traj[, i] <- rbinom(N, r1_traj[, i], theta[, 4]) 
  }
  return(list(theta = theta[, 2:4], trajs = sum_traj, r1_traj, obs_traj))
}

################################################################################
# Run code for simulated data
################################################################################

library(parallel)
set.seed(1)
system.time(
  ss_output <- ssStorvik(N=300000)
)
