################################################################################

#              Code to produce forecasts for the Ebola application             #

################################################################################

ebola_pol = function(i, state, theta, inter, dtau, sepi = stoic) {
  len <- inter / dtau + 1
  btheta <- theta[i, ]
  xmat <- matrix(0, ncol = 3, nrow = len)
  xmat[1, ] <- state[i, ]
  sumr <- matrix(0, ncol = 3, nrow = len)
  for (k in 1:(len - 1)) {
    h <- hazard(xmat[k, ], btheta)
    rnew <- rpois(3, h * dtau)
    xmat[k + 1, ] <- xmat[k, ] + sepi %*% rnew
    if (any(xmat[k + 1, ] < 0) == TRUE) {
      xmat[k + 1, which(xmat[k + 1, ] < 0)] <- 0
    }
    sumr[k + 1, ] <- sumr[k, ] + rnew
  }
  return(c(sumr[len, 2]))
}
