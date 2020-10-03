## libraries ----
library(longTransients)
library(parallel)
## set seed ----
set.seed(2020)
## parameter values ----
N <- 1e3; N_trajectories_sim <- 10
r <- 0.05; K <- 2
a <- 0.023; H <- 0.38; Q <- 5
x0 <- 0.3
sigma <- 0.02; sigma_me <- 0.05
## constants (priors) ----
constants_sim <- list(
  N = N, N_trajectories = N_trajectories_sim, 
  x0 = x0, t.step = 1, N_t = N - 1
)
## simulation parameters ----
inits_sim <- list(
  r = r, K = K, a = a, H = H, Q = Q, 
  sigma = sigma, sigma_me = sigma_me
)
## simulate data ----
sim <- simulate_model_nM(constants = constants_sim, inits = inits_sim)
x_eval <- seq(min(sim$obs_y), max(sim$obs_y), l = 2e2)
## device ----
pdf(paste0("fig/data_me_", substr(sigma_me, 3, 5), ".pdf"))
## plot data ----
plot_traj(sim$obs_y, sim$true_x)
## dev.off ----
dev.off()
## compute stable points ----
stable_pop <- uniroot(f = dpotential, interval = c(1, 2), a = inits_sim$a, r = inits_sim$r, 
                      H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K)$root
ghost_pop <- optimize(f = dpotential, interval = c(0, 1), a = inits_sim$a, r = inits_sim$r, 
                      H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K)$minimum
## ----
## simulation parallel loop ----
cl <- makeCluster(6)
clusterExport(cl, varlist = c("inits_sim", "x_eval", "sim", "N", "N_trajectories_sim", "x0", 
                              "dV"))
y_subsets <- c(lapply(1, function(x) 1:N_trajectories_sim),
               lapply(1:5, function(x) sample(1:N_trajectories_sim, 5)),
               lapply(1:5, function(x) sample(1:N_trajectories_sim, 2)),
               lapply(1:N_trajectories_sim, function(x) x))
model_compare <- parLapplyLB(cl = cl, X = y_subsets, fun = function(y_subset){
  ## libraries ----
  library(longTransients)
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
  n_iterations <- 1e5
  n_iterations_functional <- 1e5
  data <- list("y" = matrix(sim$obs_y[, y_subset], ncol = N_trajectories_fit))
  ## fit parametric ----
  fit <- fit_mcmc(data = data, constants = constants_fit, 
                  seed = seed, n_iterations = n_iterations)
  fit$sample_time
  ISE <- get_ISE(samples = fit$samples, true = inits_sim, x = x_eval)
  ## fit functional ----
  fit_functional <- fit_mcmc_functional(data = data, constants = constants_fit_functional, 
                                        seed = seed, n_iterations = n_iterations_functional)
  fit_functional$sample_time
  ISE_functional <- get_ISE_functional(samples = fit_functional$samples, true = inits_sim, 
                                       x = x_eval)
  ## device ----
  pdf(paste0("fig/posterior_trace_subset_", paste0(sort(y_subset), collapse = "_"), 
             "_me_", substr(inits_sim$sigma_me, 3, 5), ".pdf"))
  ## trace plots ----
  plot_trace(fit$samples, true = inits_sim)
  plot_trace_functional(fit_functional$samples)
  ## dev.off ----
  dev.off()
  ## device ----
  pdf(paste0("fig/posterior_potentials_subset_", paste0(sort(y_subset), collapse = "_"), 
             "_me_", substr(inits_sim$sigma_me, 3, 5), ".pdf"),
      width = 10)
  ## compare_potential curves ----
  plot_potential(samples = fit$samples, true = inits_sim, x = x_eval, 
                 obs_y = data$y, ylim_dpotential = 0.02 * c(-1, 1))
  plot_potential_functional(samples = fit_functional$samples, true = inits_sim, 
                            obs_y = data$y, x = x_eval, ylim_dpotential = 0.02 * c(-1, 1))
  ## dev.off ----
  dev.off()
  # ## compare ISE ----
  # plot(x_eval, ISE, type = "l", log = "")
  # lines(x_eval, ISE_functional, type = "l", col = 2)
  # rug(data$y)
  ## ISE diff ----
  standardized_ISE_diff <- mean(ISE - ISE_functional) / sd(ISE - ISE_functional)
  ## compare ESS ----
  ESS <- c("fit" = get_ESS(fit$samples), 
           "fit_functional" = get_ESS_functional(fit_functional$samples))
  ## return ----
  return(list("standardized_ISE_diff" = standardized_ISE_diff, "ESS" = ESS))
})
stopCluster(cl)
## ----
## save results ----
y_subset_lengths <- unlist(lapply(y_subsets, length))
ESS_by_N <- sapply(unique(y_subset_lengths), function(l){
  l_ind <- which(y_subset_lengths == l)
  rowMeans(matrix(unlist(lapply(model_compare[l_ind], function(l_i){
    l_i$ESS[c('fit.min', 'fit_functional.min')]
  })), nrow = 2))
})
colnames(ESS_by_N) <- paste(unique(y_subset_lengths), "traj")
rownames(ESS_by_N) <- c("fit", "functional")
ISE_by_N <- sapply(unique(y_subset_lengths), function(l){
  l_ind <- which(y_subset_lengths == l)
  mean(unlist(lapply(model_compare[l_ind], function(l_i)
    l_i$standardized_ISE_diff)))
})
names(ISE_by_N) <- paste(unique(y_subset_lengths), "traj")
save(ESS_by_N, ISE_by_N, file = paste0("data/ESS_ISE_me_", substr(sigma_me, 3, 5), ".RData"))
# # plot results ----
# matplot(unique(y_subset_lengths), t(ESS_by_N), type = "b", pch = 1,
#         xlab = "number of trajectories", ylab = "min(ESS)")
# legend("bottomright", lty = 1:2, col = 1:2, legend = c("parametric", "non-parametric"), bty = "n")