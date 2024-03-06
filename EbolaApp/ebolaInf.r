################################################################################

#                Functions for inference for Ebola application                 #

################################################################################

library(parallel)
library(boot)

################################################################################
#  Helper functions to generate draws from a multivariate Normal distribution  #
################################################################################

sqrtmat = function(V) {
  spec <- svd(V)
  return(spec$u%*%diag(sqrt(spec$d))%*%t(spec$u))
}

rmvn = function(m, V) {
  p <- length(m)
  z <- rnorm(p)
  return(m+sqrtmat(V)%*%z)
}


################################################################################
#     Define hazard function, stoichiometry matrix, P and condition hazard     #
################################################################################

hazard = function(state, theta) {
  # state is a length-3 vector (S,E,I)
  # rate is a length-3 vector (beta,kappa,gamma)
  c(
    theta[1]*state[1]*state[3],
    theta[2]*state[2],
    theta[3]*state[3]
  )
}

stoic <- matrix(
  c(
    -1, 0, 0,
    1, -1, 0,
    0, 1, -1
  ), ncol=3, nrow=3, byrow=TRUE
)

cond_hazard = function(state, theta, phi, obs, Tcurr, Tobs, rcurr) {
  # phi is a length-2 vector (rho,nu)
  h <- hazard(state, theta)
  mu <- phi[1]*(rcurr[2]+h[2]*(Tobs-Tcurr))
  
  out <- c(
    h[1],
    h[2] + phi[1]*h[2]*
      (phi[1]^2*h[2]*(Tobs-Tcurr)+(mu+mu^2/phi[2]))^(-1)*(obs-mu),
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

ebolaBridge = function(i, state, theta, phi, obs, Tobs, dtau, S=stoic) {
  # state, theta, phi are all be matrices with rows as in vectors above
  N <- Tobs/dtau + 1
  
  btheta <- theta[i, ]
  bphi <- phi[i, ]
  xmat <- matrix(0, ncol=3, nrow=N)
  xmat[1, ] <- state[i, ]
  
  sumr <- matrix(0, ncol=3, nrow=N)
  sumg <- matrix(0, ncol=3, nrow=N)
  sumg[1, ] <- c(xmat[1, 1]*xmat[1, 3], xmat[1, 2], xmat[1, 3])
  llr <- 0
  
  time <- 0
  count <- 1
  
  for (k in 1:(N-1)) {
    h <- hazard(xmat[k, ], btheta)
    hs <- cond_hazard(xmat[k, ], btheta, bphi, obs, (k-1)*dtau, Tobs, sumr[k, ])

    rnew <- rpois(3, hs*dtau)
    xmat[k+1, ] <- xmat[k, ] + S%*%rnew
    if (any(xmat[k+1, ] < 0) == TRUE) {
      xmat[k+1, which(xmat[k+1, ] < 0)] <- 0
      llr <- -Inf
    }
    sumr[k+1, ] <- sumr[k, ] + rnew
    sumg[k+1, ] <- sumg[k, ] +
      c(
        xmat[k+1, 1]*xmat[k+1, 3], xmat[k+1, 2], xmat[k+1, 3]
      )
    llr <- llr + sum(dpois(rnew, h*dtau, log=TRUE)) -
      sum(dpois(rnew, hs*dtau, log=TRUE))
    count <- count + 1
  }
  lwt <- llr + dnbinom(obs, mu=phi[1]*sumr[N, 2], size=bphi[2], log=TRUE)
  return(c(xmat[N, ], sumr[N, ], sumg[(N-1), ], lwt))
}

################################################################################
#                  Particle filter code for Ebola application                  #
################################################################################

ebolaStorvik = function(
  N=100000, dat=data, x0=c(44326, 15, 10), TT=1, dt=1/7,
  a=c(0.2, 5, 10, 0.85, 5), b=c(5000, 4.6, 10, 0.75, 0.1), Sm=stoic
) {
  #Initialise with prior draws
  cmat <- matrix(0, ncol=5, nrow=N)
  # cmat contains beta, omega, gamma, rho, nu
  cmat[, 1] <- rgamma(N, shape=a[1], rate=b[1])
  cmat[, 2] <- rgamma(N, shape=a[2], rate=b[2])
  cmat[, 3] <- rgamma(N, shape=a[3], rate=b[3])
  cmat[, 4] <- rnorm(N, a[4], b[4])
  cmat[, 4] <- inv.logit(cmat[, 4])
  cmat[, 5] <- rgamma(N, shape=5, rate=0.2)

  xmat <- matrix(x0, ncol=3, nrow=N, byrow=TRUE)
  
  # Set up storage of parameters and states for forecasting
  cstore <- array(dim=c(N, 5, 5))
  xstore <- array(dim=c(N, 3, 5))
  
  # Set up storage for summaries
  # 8 means followed by 8 lower quantiles followed by 8 upper quantiles
  sum_traj <- matrix(0, ncol=24, nrow=(length(dat)+1))
  sum_traj[1, ] <- c(
    apply(cbind(cmat, xmat), 2, mean),
    c(t(apply(cbind(cmat, xmat), 2, quantile, c(0.025, 0.975))))
  )
  lat_traj <- matrix(0, ncol=length(dat), nrow=N)
  obs_traj <- matrix(0, ncol=length(dat), nrow=N)
  
  sumr <- matrix(0, ncol=3, nrow=N)
  sumg <- matrix(0, ncol=3, nrow=N)
  
  # Calculate Liu-West h and a
  hjit <- sqrt(1-((3*0.99-1)/(2*0.99))^2)
  ajit <- sqrt(1-hjit^2)
  
  # Loop over all observations
  for (i in seq_along(dat)) {
    message(i)
    # Liu-West update for rho and nu
    lrho <- logit(cmat[, 4])
    lnu <- log(cmat[, 5])
    phi <- rbind(lrho, lnu)
    phibar <- apply(phi, 1, mean)
    phitild <- phi - phibar
    V = 1/N*phitild%*%t(phitild)
    for (k in 1:N) {
      m <- ajit*phi[, k] + (1-ajit)*phibar
      draw <- rmvn(m, hjit^2*V)
      lrho[k] <- draw[1]
      lnu[k] <- draw[2]
    }
    cmat[, 4] <- inv.logit(lrho)
    cmat[, 5] <- exp(lnu)
    # Propagate dynamic processes
    prop <- t(
      simplify2array(
        mclapply(
          1:N, ebolaBridge, state=xmat, theta=cmat[, 1:3], phi=cmat[, 4:5],
          obs=dat[i], Tobs=TT, dtau=dt, S=Sm, mc.cores=detectCores()-8
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
    pick <- sample(1:N, N, replace=TRUE, prob=wts)
    xmat <- xmat[pick, ]
    sumr <- sumr[pick, ]
    sumg <- sumg[pick, ]
    cmat <- cmat[pick, ]
    
    # Rejuvenate theta
    cmat[, 1] <- rgamma(N, shape=a[1]+sumr[, 1], rate=b[1]+dt*sumg[, 1])
    cmat[, 2] <- rgamma(N, shape=a[2]+sumr[, 2], rate=b[2]+dt*sumg[, 2])
    cmat[, 3] <- rgamma(N, shape=a[3]+sumr[, 3], rate=b[3]+dt*sumg[, 3])

    sum_traj[(i+1), ] <- c(
      apply(cbind(cmat, xmat), 2, mean),
      apply(cbind(cmat, xmat), 2, quantile, 0.025),
      apply(cbind(cmat, xmat), 2, quantile, 0.975)
    )
    
    lat_traj[, i] <- prop[pick, 5]
    obs_traj[, i] <- rnbinom(N, mu=cmat[, 4]*lat_traj[, i], size=cmat[, 5])
    
    # Store paras and states if at time 45, 46, 47, 48 or 49 (for forecasting)
    if (i %in% 45:49) {
      cstore[, , i-44] <- cmat
      xstore[, , i-44] <- xmat
    }
  }
  return(list(theta=cmat, trajs=sum_traj, lat_traj, obs_traj))
}
