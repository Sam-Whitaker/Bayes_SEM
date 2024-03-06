################################################################################

#              Code to produce forecasts for the Ebola application             #

################################################################################

ebolaPoL = function(i, state, theta, phi, Tobs, dtau, S=stoic){
  N <- Tobs/dtau + 1
  
  btheta <- theta[i, ]
  bphi <- phi[i, ]
  xmat <- matrix(0, ncol=3, nrow=N)
  xmat[1, ] <- state[i, ]
  
  sumr <- matrix(0, ncol=3, nrow=N)
  sumg <- matrix(0, ncol=3, nrow=N)
  sumg[1, ] <- c(xmat[1, 1]*xmat[1, 3], xmat[1, 2], xmat[1, 3])

  time <- 0
  count <- 1
  
  for (k in 1:(N-1)) {
    h <- hazard(xmat[k, ], btheta)
    rnew <- rpois(3, h*dtau)
    xmat[k+1, ] <- xmat[k, ] + S%*%rnew
    if (any(xmat[k+1, ] < 0) == TRUE) {
      xmat[k+1, which(xmat[k+1, ] < 0)] <- 0
    }
    sumr[k+1, ] <- sumr[k, ] + rnew
    sumg[k+1, ] <- sumg[k, ] + c(
      xmat[k+1, 1]*xmat[k+1, 3], xmat[k+1, 2], xmat[k+1, 3]
    )
    count <- count + 1
  }
  return(c(sumr[N, 2]))
}
