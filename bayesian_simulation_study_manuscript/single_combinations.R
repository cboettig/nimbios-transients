## libraries ----
library(nimble)
library(longTransients)
library(parallel)
## parameter values, y_subset ----
sigma_me <- 0.04
a <- 0.023
y_subset <- 1
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
seed <- 1234
n_iterations <- 1e5
n_iterations_functional <- n_iterations
data <- list("y" = matrix(sim$obs_y[, y_subset], ncol = N_trajectories_fit))
x_eval <- seq(0.2, 1.8, l = 2e2)
## fitting ----
## fit/load parametric ----
# fit <- fit_mcmc(data = data, constants = constants_fit,
#                 seed = seed, n_iterations = n_iterations)
# fit$sample_time
# save(fit, file = paste0("data/fit_parametric_", file_suffix, ".RData"))
load(file = paste0("data/fit_parametric_", file_suffix, ".RData"))
## fit/load functional ----
# fit_functional <- fit_mcmc_functional(data = data, constants = constants_fit_functional,
#                                       seed = seed, n_iterations = n_iterations_functional)
# fit_functional$sample_time
# save(fit_functional, file = paste0("data/fit_functional_", file_suffix, ".RData"))
load(file = paste0("data/fit_functional_", file_suffix, ".RData"))
## device ----
pdf(paste0("fig/posterior_trace_subset_", paste0(sort(y_subset), collapse = "_"), 
           "_me_", substr(sigma_me, 3, 5), ".pdf"))
## trace plots ----
plot_trace(fit$samples, true = inits_sim)
## dev.off ----
dev.off()
## device ----
pdf(paste0("fig/posterior_trace_functional_subset_", paste0(sort(y_subset), collapse = "_"), 
           "_me_", substr(sigma_me, 3, 5), ".pdf"))
## trace plots functional ----
plot_trace_functional(fit_functional$samples)
## dev.off ----
dev.off()
## device ----
pdf(paste0("fig/posterior_trace_subset_dpotential_", paste0(sort(y_subset), collapse = "_"), 
           "_me_", substr(sigma_me, 3, 5), ".pdf"))
## plot trace of dpotential at sample of populations [parametric] ----
x_selection <- sort(c(0.2, 0.4, 0.8, 1, 1.4, 1.6, 1.8, c(ghost_pop, stable_pop)))
dpotential_fit <- t(apply(fit$samples, 1, function(sample){
  dpotential(x = x_selection, a = sample['a'], r = sample['r'],
             H = sample['H'], Q = sample['Q'], K = sample['K'])
}))
layout(matrix(1:10, 5, 2))
cols <- scales::alpha(c(rep(1, 2), "red", rep(1, 2), "blue", rep(1, 3)), 1)
par(mar = c(4.1, 4.1, 2.1, 1.1))
for(i in 1:9){
  main_app <- ""
  if(i == 3) main_app <- "(ghost attractor)"
  if(i == 6) main_app <- "(stable attractor)"
  plot(dpotential_fit[, i], type = "l", lty = 1, col = cols[i],
       main = paste("x(t) =", round(x_selection[i], 2), main_app), 
       ylab = expression(d*mu(t)/dt), xlab = "iteration")
  abline(h = dpotential(x = x_selection[i], a = inits_sim$a, 
                        r = inits_sim$r, H = inits_sim$H, 
                        Q = inits_sim$Q, K = inits_sim$K), 
         lwd = 2, col = cols[i])
}
corrplot::corrplot(cor(dpotential_fit))
## dev.off ----
dev.off()
## device ----
pdf(paste0("fig/posterior_trace_subset_dpotential_functional_", paste0(sort(y_subset), collapse = "_"), 
           "_me_", substr(sigma_me, 3, 5), ".pdf"))
## plot trace of dpotential at sample of populations [functional] ----
dpotential_functional <- t(apply(fit_functional$samples, 1, function(sample){
  sapply(x_selection, dV, beta = sample[paste0('beta[', 1:5, ']')])
}))
layout(matrix(1:10, 5, 2))
par(mar = c(4.1, 4.1, 2.1, 1.1))
for(i in 1:9){
  main_app <- ""
  if(i == 3) main_app <- "(ghost attractor)"
  if(i == 6) main_app <- "(stable attractor)"
  plot(dpotential_functional[, i], type = "l", lty = 1, col = cols[i],
       main = paste("x(t) =", round(x_selection[i], 2), main_app), 
       ylab = expression(d*mu(t)/dt), xlab = "iteration")
  abline(h = dpotential(x = x_selection[i], a = inits_sim$a, 
                        r = inits_sim$r, H = inits_sim$H, 
                        Q = inits_sim$Q, K = inits_sim$K), 
         lwd = 2, col = cols[i])
}
corrplot::corrplot(cor(dpotential_functional))
## dev.off ----
dev.off()
## device ----
pdf(paste0("fig/posterior_potentials_subset_", paste0(sort(y_subset), collapse = "_"), 
           "_me_", substr(sigma_me, 3, 5), ".pdf"),
    width = 10)
## compare_potential curves ----
plot_potential(samples = fit$samples, true = inits_sim, x = x_eval,
               ylim_dpotential = 0.02 * c(-1, 1))
# points(c(stable_pop, ghost_pop),
#        dpotential(c(stable_pop, ghost_pop), a = inits_sim$a, r = inits_sim$r,
#                   H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K),
#        lwd = 2, col = "darkred", cex = 2)
plot_potential_functional(samples = fit_functional$samples, true = inits_sim, 
                          x = x_eval, ylim_dpotential = 0.02 * c(-1, 1))
# points(c(stable_pop, ghost_pop),
#        dpotential(c(stable_pop, ghost_pop), a = inits_sim$a, r = inits_sim$r,
#                   H = inits_sim$H, Q = inits_sim$Q, K = inits_sim$K),
#        lwd = 2, col = "darkred", cex = 2)
## dev.off ----
dev.off()
## ISE ----
ISE <- c("fit" = get_ISE(samples = fit$samples, true = inits_sim, x = x_eval),
         "fit_functional" = get_ISE_functional(samples = fit_functional$samples,
                                               true = inits_sim, x = x_eval))
## ESS ----
ESS <- c("fit" = get_ESS(fit$samples, x = x_eval),
         "fit_functional" = get_ESS_functional(fit_functional$samples, x = x_eval))
## compare ISE, ESS (Squared Error, Integrated over posterior) ----
n <- length(ISE) / 2
plot(x_eval, ISE[1:n], type = "l", log = "")
lines(x_eval, ISE[1:n + n], type = "l", col = 2)
legend("topright", lty = 1, col = 1:2, legend = c("parametric", "non-parametric"))
rug(data$y)
## ISE @ stable points; ESS median ----
x_ind <- c(which.max(x_eval > ghost_pop),
           which.max(x_eval > stable_pop))
mean_ISE_ratio <- mean(log(ISE[x_ind + n] / ISE[x_ind]))
mean_ISE_ratio
ESS_diff <- ESS['fit_functional.median'] - ESS['fit.median']
ESS_diff