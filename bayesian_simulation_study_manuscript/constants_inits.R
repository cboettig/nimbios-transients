## parameter values ----
N <- 1e3; N_trajectories <- 1
r <- 0.05; K <- 2
a <- 0.023; H <- 0.38; Q <- 5
x0 <- 0.3
sigma <- 0.02; sigma_me <- 0.05
## constants (priors) ----
constants <- list(
  N = N, N_trajectories = N_trajectories, 
  x0 = x0, t.step = 1, N_t = N - 1,
  ## prior hyperparameters
  r_shape = 2, r_rate = 10,
  K_shape = 2, K_rate = 1,
  a_shape = 2, a_rate = 10,
  H_shape = 2, H_rate = 10,
  Q_shape = 4, Q_rate = 1,
  sigma_shape = 1, sigma_rate = 100,
  sigma_me_shape = 1, sigma_me_rate = 100
)
## inits ----
inits <- list(
  r = r, K = K, a = a, H = H, Q = Q, 
  sigma = sigma, sigma_me = sigma_me
)