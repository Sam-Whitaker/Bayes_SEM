################################################################################

#              Code to produce forecasts for the Covid application             #

################################################################################

covid_pol = function(i, state, theta, phi, inter, dtau, sepi) {
  len <- inter / dtau + 1
  btheta <- theta[i, ]
  xmat <- matrix(0, ncol = 2, nrow = len)
  xmat[1, ] <- state[i, ]
  sumr <- matrix(0, ncol = 3, nrow = len)
  for (k in 1:(len-1)) {
    h <- c(
      xmat[i, 1] * xmat[i, 2] * exp(theta[1]) / 19453561,
      xmat[i, 2] * theta[2]
    )
    rnew <- rpois(3, h * dtau)
    xmat[k + 1, ] <- xmat[k, ] + sepi %*% rnew
    if (any(xmat[k + 1, ] < 0) == TRUE) {
      xmat[k + 1, which(xmat[k + 1, ] < 0)] <- 0
    }
    sumr[k + 1, ] <- sumr[k, ] + rnew
    theta[1] <- rnorm(1, theta[1], sqrt(dtau / theta[3]))
    phi[1] <- inv.logit(rnorm(1, logit(phi[1]), sqrt(dtau / theta[4])))
  }
  return(c(sumr[len, 2], phi[1]))
}