## libraries ----
library(longTransients)
library(parallel)
## set seed ----
set.seed(2020)
## parameter values [CHANGE THESE TO ADJUST POTENTIAL + MEASUREMENT ERROR] ----
N <- 1e3; N_trajectories_sim <- 10
r <- 0.05; K <- 2
a <- 0.023; H <- 0.38; Q <- 5
x0 <- 0.3
sigma <- 0.02; sigma_me <- 0.01
## constants (priors) ----
constants_sim <- list(
  N = N, N_trajectories = N_trajectories_sim, 
  x0 = x0, t.step = 1, N_t = N - 1
)
## initialize simulation ----
inits_sim <- list(
  r = r, K = K, a = a, H = H, Q = Q, 
  sigma = sigma, sigma_me = sigma_me
)
## simulate data ----
sim <- simulate_model_nM(constants = constants_sim, inits = inits_sim)
x_eval <- seq(min(sim$obs_y), max(sim$obs_y), l = 2e2)
## device ----
pdf(paste0("fig/data_single", substr(sigma_me, 3, 5), ".pdf"))
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
## fitting ----
## select subset of simulated trajectories ----
y_subset <- 1
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
           "_me_", substr(sigma_me, 3, 5), ".pdf"))
## trace plots ----
plot_trace(fit$samples, true = inits_sim)
plot_trace_functional(fit_functional$samples)
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
## compare ISE (Squared Error, Integrated over posterior) ----
plot(x_eval, ISE, type = "l", log = "")
lines(x_eval, ISE_functional, type = "l", col = 2)
legend("topright", lty = 1, col = 1:2, legend = c("parametric", "non-parametric"))
rug(data$y)
## ISE diff ----
standardized_ISE_diff <- mean(ISE - ISE_functional) / sd(ISE - ISE_functional)
## compare ESS (effective sample size) ----
ESS <- c("fit" = get_ESS(fit$samples), 
         "fit_functional" = get_ESS_functional(fit_functional$samples))
ESS
## ----
## can we recover parameters from non-parametric model? short answer: no ----
get_par <- function(beta, x_eval, loss = function(f, g) sum((f - g)^2), 
                    init = NULL, init_sd = NULL, n_reps = 1, verbose = F, 
                    ncores = 1, ...){
  args <- list(...)
  if(is.null(init) & !is.null(args$lower) & !is.null(args$upper)){
    init <- sapply(1:length(args$lower), function(i){
      runif(1, args$lower[i], args$upper[i])
    })
  }
  if(is.null(init_sd)){
    init_sd <- rep(0, length(init))
  }
  if(ncores > 1){
    require(parallel)
    cl <- makeForkCluster(ncores)
    optima <- parSapply(cl, 1:n_reps, function(i){
      optim(par = init + rnorm(length(init), sd = init_sd), fn = function(par){
        loss(dpotential(x = x_eval, a = par[1], r = par[2], H = par[3], Q = par[4], K = par[5]),
             sapply(x_eval, dV, beta = beta))
      }, ...)
    }, simplify = F)
    stopCluster(cl)
  } else {
    optima <- sapply(1:n_reps, function(i){
      optim(par = init + rnorm(length(init), sd = init_sd), fn = function(par){
        loss(dpotential(x = x_eval, a = par[1], r = par[2], H = par[3], Q = par[4], K = par[5]),
             sapply(x_eval, dV, beta = beta))
      }, ...)
    }, simplify = F)
  }
  best <- which.min(sapply(optima, function(opt) opt$value))
  if(verbose) return(optima)
  return(optima[[best]]$par)
}
get_par(beta = beta_samples[1, ], x_eval = x_eval, method = "L-BFGS-B", 
        n_reps = 10, lower = rep(1e-8, 5), upper = c(1, 1, 10, 20, 20),
        control = list(parscale = c(1e1, 1e1, 1, 1, 1)), ncores = 6)
unlist(inits_sim[c('a', 'r', 'H', 'Q', 'K')])
beta_samples <- fit_functional$samples[, grep("beta", colnames(fit_functional$samples))]
system.time({
  par_hat <- t(matrix(unlist(apply(beta_samples, 1, get_par, 
                                   x_eval = x_eval, n_reps = 30, ncores = 6,
                                   lower = rep(1e-8, 5), upper = c(1, 1, 10, 20, 20),
                                   control = list(parscale = c(1e1, 1e1, 1, 1, 1)))), nrow = 5))
})
colnames(par_hat) <- c('a', 'r', 'H', 'Q', 'K')
matplot(par_hat, type = "l", log = "y")
plot(par_hat[, 'r'], type = "l", log = "y")
matdensity <- function(mat){
  d <- apply(as.matrix(mat), 2, density)
  plot(range(sapply(1:length(d), function(i) range(d[[i]]$x))),
       range(sapply(1:length(d), function(i) range(d[[i]]$y))), 
       type = "n", xlab = "", ylab = "")
  for(i in 1:length(d)){
    lines(d[[i]], col = i)
  }
}
coda::effectiveSize(par_hat)
get_ESS(samples = par_hat, x = c(ghost_pop, stable_pop))
get_ESS(samples = par_hat, x = x_eval)
matdensity(par_hat[, c('r', 'a')])
matdensity(par_hat[, c('H')])
matdensity(par_hat[, c('K', 'Q')])
