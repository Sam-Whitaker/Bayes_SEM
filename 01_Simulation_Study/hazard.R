################################################################################
#               Define stoichiometry matrix and condition hazard               #
################################################################################

stoic <- matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE)

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