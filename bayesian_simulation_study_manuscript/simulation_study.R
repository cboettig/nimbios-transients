## libraries ----
library(parallel)
library(longTransients)
## parameter values ----
n_iterations <- 1e5
N <- 1e3; N_trajectories_sim <- 10
r <- 0.05; K <- 2
H <- 0.38; Q <- 5
sigma_mes <- c(0.005, 0.01, 0.02, 0.04, 0.08)
as <- c(0.023)
## define combinations of parameters ----
y_subsets <- c(lapply(1, function(x) 1:N_trajectories_sim),
               lapply(1:5, function(x) sample(1:N_trajectories_sim, 5)),
               lapply(1:5, function(x) sample(1:N_trajectories_sim, 2)),
               lapply(1:N_trajectories_sim, function(x) x))
combos <- expand.grid("y_subset" = y_subsets, "sigma_me" = sigma_mes, "a" = as)
x_eval <- seq(0.2, 1.8, l = 2e2)
## simulation parallel loop ----
cl <- makeCluster(4)
clusterExport(cl, varlist = c("N", "N_trajectories_sim", "n_iterations", "x_eval"))
model_compare <- parApply(cl = cl, X = combos, MARGIN = 1, FUN = function(combo){
# model_compare <- apply(X = combos, MARGIN = 1, FUN = function(combo){
  ## libraries ----
  library(longTransients)
  ## y_subset, sigma_me, a ----
  y_subset <- combo$y_subset
  sigma_me <- combo$sigma_me
  a <- combo$a
  ## load data ----
  y_subset_ind <- which(names(combo) == "y_subset")
  nice_value_names <- sapply(combo, function(x) gsub(".", "_", x, fixed = T))
  nice_value_names[y_subset_ind] <- paste0(y_subset, collapse = "_")
  file_suffix <- paste(names(combo), nice_value_names, collapse = "_", sep = "_")
  message("file_suffix:", file_suffix)
  load(paste0("data/sim_", paste(names(combo)[-y_subset_ind], nice_value_names[-y_subset_ind],
                                 collapse = "_", sep = "_"), ".Rdata"))
  x0 <- constants_sim$x0
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
  ## ISE ----
  tryCatch(expr = {
    ISE <- NULL
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
    ESS <- NULL
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
save(model_compare, combos, file = paste0("data/ESS_ISE_a_", 
                                          gsub(".", "_", as, fixed = T), ".RData"))
## load results ----
load(paste0("data/ESS_ISE_a_", gsub(".", "_", as, fixed = T), ".RData"))
## aggregate results ----
problem_combos <- which(unlist(lapply(model_compare, function(m) sum(is.na(unlist(m))))) > 0)
# combos[problem_combos, ]
# model_compare[[problem_combos]]$ESS
mean_ISE_ratio <- sapply(1:length(model_compare), function(i){
  model <- model_compare[[i]]
  n <- length(model$ISE) / 2
  ## need to sort out a for general combinations [20201012HRS]
  stable_pop <- uniroot(f = dpotential, interval = c(1, 2), a = combos[i, 'a'], r = r,
                        H = H, Q = Q, K = K)$root
  ghost_pop <- optimize(f = dpotential, interval = c(0, 1), a = combos[i, 'a'], r = r,
                        H = H, Q = Q, K = K)$minimum
  x_ind <- c(which.max(x_eval > ghost_pop),
             which.max(x_eval > stable_pop))
  mean(log(model$ISE[x_ind + n] / model$ISE[x_ind]))
})
ESS_diff <- unlist(sapply(model_compare, function(model){
  n <- length(model$ISE) / 2
  out <- model$ESS['fit_functional.median'] - model$ESS['fit.median']
  if(length(out) == 0) return(NA)
  return(out)
}))
summary_df <- data.frame(n_traj = apply(combos, 1, function(x) length(x$y_subset)[[1]]),
                         mean_ISE_ratio = mean_ISE_ratio,
                         ESS_diff = ESS_diff,
                         me = apply(combos, 1, function(x) x$sigma_me))
# agg_summary_df <- aggregate(summary_df[, 2:3], 
#                             list(n_traj = summary_df$n_traj, 
#                                  me = summary_df$me), median, na.rm = T)
agg_summary_df <- aggregate(summary_df[, 2:3], 
                            list(n_traj = summary_df$n_traj, 
                                 me = summary_df$me), quantile, na.rm = T, 
                            prob = c(0.25, 0.5, 0.75))
## device ----
pdf(paste0("fig/ESS_ISE_summary_plots_a_", gsub(".", "_", as, fixed = T), ".pdf"))
## plot ----
layout(matrix(1:2, 1, 2))
traj_colors <- viridisLite::magma(4, end = 0.9)
matplot(unique(agg_summary_df$me), 
        t(matrix(agg_summary_df[, 'ESS_diff'][, 2], nrow = 4, ncol = 5)), type = "l",
        ylim = c(0, max(agg_summary_df[, 'ESS_diff'])), col = traj_colors, lwd = 3,
        ylab = "diff. in ESS (non-parametric - parametric)",
        xlab = "measurement error variance", xaxt = "n", ylog = T)
for(n_traj in unique(agg_summary_df$n_traj)){
  polygon(c(unique(agg_summary_df$me), rev(unique(agg_summary_df$me))), 
          c(agg_summary_df$ESS_diff[agg_summary_df$n_traj == n_traj, 1],
            rev(agg_summary_df$ESS_diff[agg_summary_df$n_traj == n_traj, 3])),
          col = scales::alpha(traj_colors[which(n_traj == unique(agg_summary_df$n_traj))], 0.3), 
          border = NA)
}
axis(1, unique(agg_summary_df$me))
legend("topright", lty = 1:4, col = traj_colors, lwd = 3, bty = "n",
       legend = unique(agg_summary_df$n_traj))

matplot(unique(agg_summary_df$me), 
        t(matrix(agg_summary_df[, 'mean_ISE_ratio'][, 2], nrow = 4, ncol = 5)), type = "l",
        ylim = 1.55 * c(-1, 1),
        col = traj_colors, lwd = 3,
        ylab = "ratio of MISE (non-parametric / parametric)",
        xlab = "measurement error variance", xaxt = "n", ylog = T)
for(n_traj in unique(agg_summary_df$n_traj)){
  polygon(c(unique(agg_summary_df$me), rev(unique(agg_summary_df$me))), 
          c(agg_summary_df$mean_ISE_ratio[agg_summary_df$n_traj == n_traj, 1],
            rev(agg_summary_df$mean_ISE_ratio[agg_summary_df$n_traj == n_traj, 3])),
          col = scales::alpha(traj_colors[which(n_traj == unique(agg_summary_df$n_traj))], 0.3), 
          border = NA)
}
axis(1, unique(agg_summary_df$me))
# for(n_traj in unique(agg_summary_df$n_traj)){
#   mat <- t(matrix(summary_df[summary_df$n_traj == n_traj, 
#                              'mean_ISE_ratio'], ncol = 5))
#   matplot(unique(agg_summary_df$me), mat, 
#           add = T, type = "l", lty = which(n_traj == unique(agg_summary_df$n_traj)),
#           col = scales::alpha(traj_colors[which(n_traj == unique(agg_summary_df$n_traj))], 
#                               alpha = 0.2))
# }
abline(h = 0)
legend("topright", lty = 1:4, col = traj_colors, lwd = 3, bty = "n",
       legend = unique(agg_summary_df$n_traj))

# image(unique(agg_summary_df$n_traj), unique(agg_summary_df$me), 
#       zlim = c(min(agg_summary_df[, 'ESS_diff']), 0),
#       matrix(agg_summary_df[, 'ESS_diff'], nrow = 4, ncol = 5),
#       axes = F, xlab = "number of trajectories", ylab = "measurement error",
#       main = "difference in ESS", col = viridisLite::viridis(1e2, option = "A"))
# axis(1, unique(agg_summary_df$n_traj), lty = 0)
# axis(2, unique(agg_summary_df$me), lty = 0)
# image(unique(agg_summary_df$n_traj), unique(agg_summary_df$me), 
#       matrix(agg_summary_df[, 'mean_ISE_ratio'], nrow = 4, ncol = 5),
#       axes = F, xlab = "number of trajectories", ylab = "", 
#       main = "ratio of MISE", col = viridisLite::viridis(1e2))
# axis(1, unique(agg_summary_df$n_traj), lty = 0)
## dev.off ----
dev.off()