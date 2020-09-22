## load samples ----
source("code/scharf_nimble/model.R")
seed <- 270; N_trajectories <- 10
load(paste0("data/scharf_nimble/samples_functional_", N_trajectories, "_", seed, ".RData"))
samples <- fit$samples
## device ----
pdf(file = paste0("fig/trace_plots_", N_trajectories, "_", seed, ".pdf"))
# ## trace plots (log) ----
# layout(matrix(c(1:7, 7), 4, 2))
# par(mar = c(2, 2, 4, 2))
# plot(exp(samples[, 'log_r']), type = "l", main = "r", ylab = "", log = "y")
# abline(h = r)
# plot(exp(samples[, 'log_K']), type = "l", main = "K", ylab = "", log = "y")
# abline(h = K)
# plot(exp(samples[, 'log_a']), type = "l", main = "a", ylab = "",
#      ylim = range(exp(samples[, 'log_a']), a), log = "y")
# abline(h = a)
# plot(exp(samples[, 'log_H']), type = "l", main = "H", ylab = "", log = "y")
# abline(h = H)
# plot(exp(samples[, 'log_Q']), type = "l", main = "Q", ylab = "", log = "y")
# abline(h = Q)
# plot(exp(samples[, 'log_sigma']), type = "l", main = "sigma", ylab = "", log = "y")
# abline(h = sigma)
# corrplot::corrplot(cor(samples[, c("log_r", "log_K", "log_a", "log_H", "log_Q", "log_sigma")]))
## trace plots ----
layout(matrix(c(1:3, 3), 2, 2), widths = c(1, 2))
variables <- c(paste0("beta[", 1:5, "]"), "sigma", "sigma_me")
matplot(samples[, paste0("beta[", 1:5, "]")], type = "l")
matplot(samples[, c("sigma", "sigma_me")], type = "l", lty = 1)
legend("bottomright", col = 1:2, lty = 1, bty = "n", 
       legend = c(expression(sigma), expression(sigma[me])))
corrplot::corrplot(cor(samples[, variables]))
## prior/post densities ----
layout(matrix(c(1:8), 4, 2))
par(mar = c(2, 2, 4, 2))
for(var in variables[1:5]){
  plot(density(samples[, var]), main = var,
       xlim = range(0, samples[, var]))
  xx <- seq(par()$usr[1], par()$usr[2], l = 2e2)
  lines(xx, dnorm(xx, 0, 1), col = "gray", lty = 2)
}
for(var in variables[6:7]){
  plot(density(samples[, var], from = 0), main = var,
       xlim = range(0, get(var), samples[, var]))
  xx <- seq(par()$usr[1], par()$usr[2], l = 2e2)
  lines(xx, dgamma(xx, as.numeric(constants[paste0(var, "_shape")]),
                   as.numeric(constants[paste0(var, "_rate")])),
        col = "gray", lty = 2)
  points(get(var), par()$usr[3], xpd = T, pch = 8)
}
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("fig/posterior_potentials_", N_trajectories, "_", seed, ".pdf"), width = 10)
## potential curves ----
growth <- function(x, r, K){x * r * (1 - x / K)}
consumption <- function(x, a, H, Q){a * x^Q / (x^Q + H^Q)}
x <- seq(0, 1.5, l=1e2)
dpotential <- function(x = seq(0, 2, l = 1e2), a, r, H, Q, K){
  growth(x = x, r = r, K = K) - 
    consumption(x = x, a = a, H = H, Q = Q)
}
potential <- function(x = seq(0, 2, l = 1e2), a, r, H, Q, K){
  cumsum(-dpotential(x = x, r = r, K = K, a = a, H = H, Q = Q)[-1] * diff(x))
}
dpotential_curves <- apply(samples[, 1:degree], 1, function(row){
  sapply(x, get_dV, beta = row)
})
dpotential_curves_quantiles <- apply(dpotential_curves, 1, quantile, probs = c(0.125, 0.5, 0.875))
potential_curves <- apply(-dpotential_curves[-1, ] * diff(x), 2, cumsum)
potential_curves_quantiles <- apply(potential_curves, 1, quantile, probs = c(0.125, 0.5, 0.875))
subset <- sample(1:nrow(samples), min(400, nrow(samples)))
layout(matrix(1:2, 1, 2))
transparency <- 2e-2
matplot(x[-1], t(potential_curves_quantiles), type = "l", 
        lty = c(2, 1, 2), lwd = c(1, 2, 1), 
        ylim = c(-0.015, 0.01), main = "potential function", 
        ylab = "", xlab = "population", col = "black")
matplot(x[-1], potential_curves[, subset], type = "l", lty = 1, 
        col = scales::alpha("black", transparency), lwd = 2, add = T)
lines(x[-1], potential(x = x, a = a, r = r, H = H, Q = Q, K = K), lwd = 2, col = "darkred")
rug(data$y)
matplot(x, -t(dpotential_curves_quantiles), type = "l", 
        lty = c(2, 1, 2), lwd = c(1, 2, 1), 
        ylim = c(-0.05, 0.05),
        main = "derivative of potential function", 
        ylab = "", xlab = "population", col = "black")
matplot(x, -dpotential_curves[, subset], type = "l", lty = 1,
        col = scales::alpha("black", transparency), lwd = 2, add = T)
lines(x, -dpotential(x = x, a = a, r = r, H = H, Q = Q, K = K), lwd = 2, col = "darkred")
abline(h = 0, lwd = 2, lty = 3)
rug(data$y)
## ESS ----
library(coda)
effectiveSize(as.mcmc(t(potential_curves[seq(1, nrow(potential_curves), l = 1e2), ])))
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("fig/posterior_population_trajectories_", 
                  N_trajectories, "_", seed, ".pdf"), width = 10)
## true populations ----
x_chains <- array(samples[, grep('x', colnames(samples))], 
                  dim = c(nrow(samples), N, N_trajectories))
samples_sub <- sort(sample(1:nrow(samples), min(nrow(samples), 5e2)))
layout(1)
matplot(seq(1, N, l = N/constants$t.step), t(x_chains[samples_sub, , 1]),
        type = "n", ylab = "x", xlab = "time", ylim = range(x_chains))
for(n in 1:N_trajectories){
  matplot(seq(1, N, l = N/constants$t.step), t(x_chains[samples_sub, , n]),
          type = "l", add = T, col = scales::alpha("black", 0.01), lty = 1)
}
matplot(seq(1, N, l = N/constants$t.step), cmodel$x, type = "l",
      col = scales::alpha("darkred", 1), lty = 1, lwd = 2, add = T)
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("fig/posterior_deterministic_core_", N_trajectories, "_", seed, ".pdf"), width = 10)
## deterministic core ----
core_curve <- function(a, r, H, Q, K, x0 = 0.3, t = (1:1e4) * 0.5, t.step = 0.5){
  x <- rep(x0, length(t))
  for(t_i in 1:(length(t) - 1)){
    x[t_i + 1] <- x[t_i] + t.step * (x[t_i] * r * (1 - x[t_i] / K)  - a * x[t_i] ^ Q / (x[t_i] ^ Q + H ^ Q))
  }
  return(cbind(t, x))
}
subset <- sample(1:nrow(samples), min(400, nrow(samples)))
deterministic_core <- apply(samples[subset, ], 1, function(row){
  core_curve(a = row['a'], r = row['r'], 
             H = row['H'], Q = row['Q'], K = row['K'])[, 2]
})
layout(1)
matplot((1:1e4)/2, deterministic_core, type = "l", lty = 1, 
        col = scales::alpha("black", 1e-1), lwd = 2, 
        ylim = c(0.2, min(max(deterministic_core, 1.5, na.rm = T), 5)),
        main = "potential function", ylab = "population", xlab = "time")
lines(core_curve(a, r, H, Q, K), lwd = 3, lty = 2, col = "darkgreen")
legend("topleft", lwd = 2, lty = c(1, 2), col = c(scales::alpha("black", 1e-1), "darkgreen"), 
       legend = c("based on posterior draw", "true params"))
## dev.off ----
dev.off()
# ## transition time ----
# deterministic_core_ <- apply(samples[subset, ], 1, function(row){
#   core_curve(t = (1:2e4) * 0.5, a = exp(row['log_a']), r = exp(row['log_r']), 
#              H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))[, 2]
# })
# median(apply(deterministic_core_, 2, function(x){
#   if(max(x) < 1) 
#     1e4
#   else 
#     which.max(x > 1)
#   }))