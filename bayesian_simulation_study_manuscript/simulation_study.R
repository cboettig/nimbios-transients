## libraries ----
library(parallel)
## parameter values ----
n_iterations <- 1e5
N <- 1e3; N_trajectories_sim <- 10
sigma_mes <- c(0.005, 0.01, 0.02, 0.04, 0.08)
as <- c(0.023)
## define combinations of parameters ----
y_subsets <- c(lapply(1, function(x) 1:N_trajectories_sim),
               lapply(1:5, function(x) sample(1:N_trajectories_sim, 5)),
               lapply(1:5, function(x) sample(1:N_trajectories_sim, 2)),
               lapply(1:N_trajectories_sim, function(x) x))
combos <- expand.grid("y_subset" = y_subsets, "sigma_me" = sigma_mes, "a" = as)
## simulation parallel loop ----
cl <- makeCluster(6)
clusterExport(cl, varlist = c("N", "N_trajectories_sim", "n_iterations"))
model_compare <- parApply(cl = cl, X = combos, MARGIN = 1, FUN = function(combo){
# model_compare <- apply(X = combos, MARGIN = 1, FUN = function(combo){
  ## libraries ----
  library(longTransients)
  ## y_subset, sigma_me, a ----
  y_subset <- combo$y_subset[[1]]
  sigma_me <- combo$sigma_me
  a <- combo$a
  ## load data ----
  y_subset_ind <- which(names(combo) == "y_subset")
  nice_value_names <- sapply(combo, function(x) gsub(".", "_", x, fixed = T))
  nice_value_names[y_subset_ind] <- paste0(y_subset, collapse = "_")
  file_suffix <- paste(names(combo), nice_value_names, collapse = "_", sep = "_")
  load(paste0("data/sim_", paste(names(combo)[-y_subset_ind], nice_value_names[-y_subset_ind], 
                                 collapse = "_", sep = "_"), ".Rdata"))
  x0 <- constants_sim$x0
  ## x_eval + compute stable points ----
  x_eval <- seq(0, 2, l = 2e2)
  stable_pop <- uniroot(f = dpotential, interval = c(1, 2), a = inits_sim$a, r = inits_sim$r, 
                        H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K)$root
  ghost_pop <- optimize(f = dpotential, interval = c(0, 1), a = inits_sim$a, r = inits_sim$r, 
                        H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K)$minimum
  ## constants (priors) ----
  N_trajectories_fit <- length(y_subset)
  constants_fit <- list(
    N = N, N_trajectories = N_trajectories_fit, 
    x0 = x0, t.step = 1, N_t = N - 1,
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
    N = N, N_trajectories = N_trajectories_fit, 
    x0 = x0, t.step = 1, N_t = N - 1,
    degree = degree,
    beta_mean = rep(0, degree),
    beta_prec = 1e-1 * diag(degree),
    sigma_shape = 1, sigma_rate = 10,
    sigma_me_shape = 1, sigma_me_rate = 10
  )
  ## seed, data, x_eval, n_iterations ----
  seed <- 1234
  n_iterations <- n_iterations
  n_iterations_functional <- n_iterations
  data <- list("y" = matrix(sim$obs_y[, y_subset], ncol = N_trajectories_fit))
  ## fit parametric ----
  fit <- fit_mcmc(data = data, constants = constants_fit, 
                  seed = seed, n_iterations = n_iterations)
  fit$sample_time
  ## fit functional ----
  fit_functional <- fit_mcmc_functional(data = data, constants = constants_fit_functional, 
                                        seed = seed, n_iterations = n_iterations_functional)
  fit_functional$sample_time
  ## device ----
  pdf(paste0("fig/posterior_trace_", file_suffix, ".pdf"))
  ## trace plots ----
  plot_trace(fit$samples, true = inits_sim)
  plot_trace_functional(fit_functional$samples)
  ## dev.off ----
  dev.off()
  ## device ----
  pdf(paste0("fig/posterior_potentials_", file_suffix, ".pdf"),
      width = 10)
  ## compare_potential curves ----
  tryCatch(
    expr = plot_potential(samples = fit$samples, true = inits_sim, x = x_eval, 
                          obs_y = data$y, ylim_dpotential = 0.02 * c(-1, 1)),
    error = function(err){
      message("error in plot_potential() for ", file_suffix) 
      return(NULL)
    })
  tryCatch(
    expr = plot_potential_functional(samples = fit_functional$samples, true = inits_sim, 
                                     obs_y = data$y, x = x_eval, 
                                     ylim_dpotential = 0.02 * c(-1, 1)),
    error = function(err){
      message("error in plot_potential_functional() for ", file_suffix) 
      return(NULL)
    })
  ## dev.off ----
  dev.off()
  ## ISE diff ----
  tryCatch(expr = {
    ISE <- c("fit" = get_ISE(samples = fit$samples, true = inits_sim, x = x_eval),
             "fit_functional" = get_ISE_functional(samples = fit_functional$samples, 
                                         true = inits_sim, x = x_eval))
  },
  error = function(err){
    message("error in get_ISE() or get_ISE_functional() for ", file_suffix) 
    return(NULL)
  })
  ## ESS ----
  tryCatch(expr = {
    ESS <- c("fit" = get_ESS(fit$samples, x = x_eval), 
             "fit_functional" = get_ESS_functional(fit_functional$samples, x = x_eval))
    },
    error = function(err){
      message("error in get_ESS() or get_ESS_functional() for ", file_suffix) 
      return(NULL)
    })
  ## return ----
  return(list("ISE" = ISE, "ESS" = ESS))
})
stopCluster(cl)
## ----
## save results ----
save(model_compare, combos, file = paste0("data/ESS_ISE.RData"))
# y_subset_lengths <- unlist(lapply(y_subsets, length))
# ESS_by_N <- sapply(unique(y_subset_lengths), function(l){
#   l_ind <- which(y_subset_lengths == l)
#   rowMeans(matrix(unlist(lapply(model_compare[l_ind], function(l_i){
#     l_i$ESS[c('fit.min', 'fit_functional.min')]
#   })), nrow = 2))
# })
# colnames(ESS_by_N) <- paste(unique(y_subset_lengths), "traj")
# rownames(ESS_by_N) <- c("fit", "functional")
# ISE_by_N <- sapply(unique(y_subset_lengths), function(l){
#   l_ind <- which(y_subset_lengths == l)
#   mean(unlist(lapply(model_compare[l_ind], function(l_i)
#     l_i$standardized_ISE_diff)))
# })
# names(ISE_by_N) <- paste(unique(y_subset_lengths), "traj")
# save(ESS_by_N, ISE_by_N, file = paste0("data/ESS_ISE_me_", substr(sigma_me, 3, 5), ".RData"))
# # plot results ----
# matplot(unique(y_subset_lengths), t(ESS_by_N), type = "b", pch = 1,
#         xlab = "number of trajectories", ylab = "min(ESS)")
# legend("bottomright", lty = 1:2, col = 1:2, legend = c("parametric", "non-parametric"), bty = "n")