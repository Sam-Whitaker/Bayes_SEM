################################################################################
#      Define hazard function, stoichiometry matrix, and condition hazard      #
################################################################################

hazard = function(state, theta) {
  # Input:
  #  - state: a length-3 vector (S,E,I)
  #  - theta: a length-3 vector (beta,kappa,gamma)
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
  # Input
  #  - state: a length-3 vector (S,E,I)
  #  - theta: a length-3 vector (beta,kappa,gamma)
  #  - phi: a length-2 vector (rho,nu)
  #  - obs: the next observation
  #  - curr: the current simulated time
  #  - inter: the inter-observation time
  #  - rcurr: the current simulated incidence
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