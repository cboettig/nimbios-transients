## parameter values ----
N <- 1e3; N_trajectories_sim <- 10
r <- 0.05; K <- 2
H <- 0.38; Q <- 5
x0 <- 0.3
sigma <- 0.02
sigma_mes <- c(0.002, 0.01, 0.02, 0.04, 0.2)
as <- c(0.023)
parameter_combos <- expand.grid("sigma_me" = sigma_mes, "a" = as)
## simulate data loop ----
for(i in 1:nrow(parameter_combos)){
  ## libraries ----
  library(longTransients)
  ## set seed ----
  set.seed(2020)
  ## constants ----
  constants_sim <- list(
    N = N, N_trajectories = N_trajectories_sim, 
    x0 = x0, t.step = 1, N_t = N - 1
  )
  ## simulation parameters ----
  inits_sim <- list(
    r = r, K = K, a = parameter_combos[i, 'a'], H = H, Q = Q, 
    sigma = sigma, sigma_me = parameter_combos[i, 'sigma_me']
  )
  ## simulate_model_nM ----
  sim <- simulate_model_nM(constants = constants_sim, inits = inits_sim)
  ## save sim ----
  nice_value_names <- sapply(parameter_combos[i, ], function(x) gsub(".", "_", x, fixed = T))
  save(sim, inits_sim, constants_sim,
       file = paste0("data/sim_", paste(names(parameter_combos), nice_value_names, 
                                        collapse = "_", sep = "_"), ".RData"))
  ## device ----
  pdf(paste0("fig/data_", paste(names(parameter_combos), nice_value_names, 
                                collapse = "_", sep = "_"), ".pdf"))
  ## plot data ----
  plot_traj(sim$obs_y, sim$true_x)
  ## dev.off ----
  dev.off()
}