################################################################################
#              Define stoichiometry matrix, and condition hazard               #
################################################################################

stoic <- matrix(
  c(
    -1, 0,
    1, -1
  ), nrow = 2, byrow = TRUE
)

cond_hazard = function(state, theta, phi, obs, curr, inter, rcurr) {
  # Input
  #  - state: a length-2 vector (S,I)
  #  - theta: a length-4 vector (log(beta_t),kappa,lambda_beta,lambda_rho)
  #  - phi: a length-2 vector (rho_t,nu)
  #  - obs: the next observation
  #  - curr: the current simulated time
  #  - inter: the inter-observation time
  #  - rcurr: the current simulated incidence
  h <- c(state[1] * state[2] * exp(theta[1]) / 19453561, state[2] * theta[2])
  mu <- phi[1] * (rcurr[1] + h[1] * (inter - curr))
  out <- c(
    h[1] + phi[1] * h[1] *
      (phi[1]^2 * h[1] * (inter - curr) + mu + (phi[2] * mu)^2)^(-1) *
      (obs - mu),
    h[2]
  )
  out[out <= 0] <- 1e-6
  return(out)
}