## parameter values ----
N <- 1e3; N_trajectories <- 1e0
r <- 0.05; K <- 2
a <- 0.023; H <- 0.38; Q <- 5
x0 <- 0.3
sigma <- 0.02; sigma_me <- 0.01
## constants (priors) ----
constants <- list(
  N = N, N_trajectories = N_trajectories, 
  x0 = x0, t.step = 1, N_t = N - 1,
  ## prior hyperparameters
  # mu_r = log(r), sd_r = 1,
  # mu_K = log(K), sd_K = 1,
  # mu_a = log(a), sd_a = 1,
  # mu_H = log(H), sd_H = 1,
  # mu_Q = log(Q), sd_Q = 1,
  # mu_sigma = log(sigma), sd_sigma = 1,
  # mu_sigma_me = log(sigma_me), sd_sigma_me = 1,
  r_shape = 2, r_rate = 10,
  K_shape = 2, K_rate = 1,
  a_shape = 2, a_rate = 10,
  H_shape = 2, H_rate = 1,
  Q_shape = 2, Q_rate = 1/2,
  sigma_shape = 1, sigma_rate = 10,
  sigma_me_shape = 1, sigma_me_rate = 10
)
## inits ----
inits <- list(
  # log_r = log(r),
  # log_K = log(K),
  # log_a = log(a),
  # log_H = log(H),
  # log_Q = log(Q), log_sigma = log(sigma),
  # log_sigma_me = log(sigma_me)
  r = r, K = K, a = a, H = H, Q = Q, 
  sigma = sigma, sigma_me = sigma_me
)