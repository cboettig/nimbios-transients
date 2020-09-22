degree <- 5
## constants (priors) ----
constants_functional <- list(
  N = N, N_trajectories = N_trajectories, 
  x0 = x0, t.step = 1, N_t = N - 1,
  degree = degree,
  beta_mean = rep(0, degree),
  beta_prec = 1e-1 * diag(degree),
  sigma_shape = 1, sigma_rate = 10,
  sigma_me_shape = 1, sigma_me_rate = 10
)