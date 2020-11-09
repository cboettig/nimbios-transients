## libraries ----
library(longTransients)
library(parallel)
## parameter values, y_subset ----
N <- 1e5; N_trajectories_sim <- 10
r <- 0.05; K <- 2
H <- 0.38; Q <- 5
sigma_me <- 0.01
a <- 0.023
n_trajs <- c(1, 2, 5, 10)
## fit each number of traj once ----
for(n_traj in n_trajs){
  y_subset <- 1:n_traj
  combo <- list(sigma_me = sigma_me, a = a, y_subset = y_subset)
  ## load data ----
  y_subset_ind <- which(names(combo) == "y_subset")
  nice_value_names <- sapply(combo, function(x) gsub(".", "_", x, fixed = T))
  nice_value_names[y_subset_ind] <- paste0(y_subset, collapse = "_")
  file_suffix <- paste(names(combo), nice_value_names, collapse = "_", sep = "_")
  message("file_suffix:", file_suffix)
  load(paste0("data/sim_", paste(names(combo)[-y_subset_ind], nice_value_names[-y_subset_ind],
                                 collapse = "_", sep = "_"), ".Rdata"))
  x0 <- constants_sim$x0
  ## compute stable points ----
  stable_pop <- uniroot(f = dpotential, interval = c(1, 2), a = inits_sim$a, 
                        r = inits_sim$r, H = inits_sim$H, Q = inits_sim$Q, 
                        K = inits_sim$K)$root
  ghost_pop <- optimize(f = dpotential, interval = c(0, 1), a = inits_sim$a, 
                        r = inits_sim$r, H = inits_sim$H, Q = inits_sim$Q, 
                        K = inits_sim$K)$minimum
  ## ----
  ## constants (priors) ----
  N_trajectories_fit <- length(y_subset)
  constants_fit <- list(
    N = constants_sim$N, N_trajectories = N_trajectories_fit,
    x0 = x0, t.step = 1, N_t = constants_sim$N - 1,
    ## prior hyperparameters
    r_shape = 2, r_rate = 10,
    K_shape = 1, K_rate = 1e-1,
    a_shape = 2, a_rate = 10,
    H_shape = 2, H_rate = 1,
    Q_shape = 1, Q_rate = 1e-1,
    sigma_shape = 1, sigma_rate = 10,
    sigma_me_shape = 1, sigma_me_rate = 10
  )
  ## constants for functional model (priors) ----
  degree <- 5
  constants_fit_functional <- list(
    N = constants_sim$N, N_trajectories = N_trajectories_fit,
    x0 = x0, t.step = 1, N_t = constants_sim$N - 1,
    degree = degree,
    beta_mean = rep(0, degree),
    beta_prec = 1e-1 * diag(degree),
    sigma_shape = 1, sigma_rate = 10,
    sigma_me_shape = 1, sigma_me_rate = 10
  )
  ## seed, data, x_eval, n_iterations ----
  seed <- 2020
  n_iterations <- 1e5
  n_iterations_functional <- n_iterations
  data <- list("y" = matrix(sim$obs_y[, y_subset], ncol = N_trajectories_fit))
  x_eval <- seq(0.2, 1.8, l = 2e2)
  ## fitting ----
  ## fit parametric ----
  fit <- fit_mcmc(data = data, constants = constants_fit,
                  seed = seed, n_iterations = n_iterations)
  fit$sample_time
  save(inits_sim, constants_fit, fit, file = paste0("data/fit_parametric_", file_suffix, ".RData"))
  # ## fit functional ----
  # fit_functional <- fit_mcmc_functional(data = data, constants = constants_fit_functional,
  #                                       seed = seed, n_iterations = n_iterations_functional)
  # fit_functional$sample_time
  # save(inits_sim, constants_fit_functional, fit_functional, 
  #      file = paste0("data/fit_functional_", file_suffix, ".RData"))
  # ## device ----
  # pdf(paste0("fig/posterior_potentials_subset_", paste0(sort(y_subset), collapse = "_"), 
  #            "_me_", substr(sigma_me, 3, 5), ".pdf"),
  #     width = 10)
  # ## compare_potential curves ----
  # plot_potential(samples = fit$samples, true = inits_sim, x = x_eval,
  #                ylim_dpotential = 0.02 * c(-1, 1))
  # # points(c(stable_pop, ghost_pop),
  # #        dpotential(c(stable_pop, ghost_pop), a = inits_sim$a, r = inits_sim$r,
  # #                   H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K),
  # #        lwd = 2, col = "darkred", cex = 2)
  # plot_potential_functional(samples = fit_functional$samples, true = inits_sim, 
  #                           x = x_eval, ylim_dpotential = 0.02 * c(-1, 1))
  # # points(c(stable_pop, ghost_pop),
  # #        dpotential(c(stable_pop, ghost_pop), a = inits_sim$a, r = inits_sim$r,
  # #                   H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K),
  # #        lwd = 2, col = "darkred", cex = 2)
  # ## dev.off ----
  # dev.off()
}