# ## load samples ----
# source("model.R")
# seed <- 270; N_trajectories <- 10
# load(paste0("../../data/scharf_nimble/samples_", N_trajectories, "_", seed, ".RData"))
## device ----
pdf(file = paste0("../../figs/scharf_nimble/trace_plots_", N_trajectories, "_", seed, ".pdf"))
## trace plots ----
layout(matrix(c(1:7, 7), 4, 2))
par(mar = c(2, 2, 4, 2))
plot(exp(samples[, 'log_r']), type = "l", main = "r", ylab = "", log = "y")
abline(h = r)
plot(exp(samples[, 'log_K']), type = "l", main = "K", ylab = "", log = "y")
abline(h = K)
plot(exp(samples[, 'log_a']), type = "l", main = "a", ylab = "",
     ylim = range(exp(samples[, 'log_a']), a), log = "y")
abline(h = a)
plot(exp(samples[, 'log_H']), type = "l", main = "H", ylab = "", log = "y")
abline(h = H)
plot(exp(samples[, 'log_Q']), type = "l", main = "Q", ylab = "", log = "y")
abline(h = Q)
plot(exp(samples[, 'log_sigma']), type = "l", main = "sigma", ylab = "", log = "y")
abline(h = sigma)
corrplot::corrplot(cor(samples[, c("log_r", "log_K", "log_a", "log_H", "log_Q", "log_sigma")]))
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../../figs/scharf_nimble/posterior_potentials_", N_trajectories, "_", seed, ".pdf"), width = 10)
## potential curves ----
growth <- function(x, r, K){x * r * (1 - x / K)}
consumption <- function(x, a, H, Q){a * x^Q / (x^Q + H^Q)}
x <- seq(0, 2, l=1e2)
potential <- function(x = seq(0, 2, l = 1e2), a, r, H, Q, K){
  - cumsum(growth(x = x, r = r, K = K) - 
             consumption(x = x, a = a, H = H, Q = Q))
}
potential_curves <- apply(samples, 1, function(row){
  potential(x = x, a = exp(row['log_a']), r = exp(row['log_r']), 
            H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))
})
dpotential_curves <- apply(samples, 1, function(row){
  p <- potential(x = x, a = exp(row['log_a']), r = exp(row['log_r']), 
            H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))
  diff(p) / diff(x)
})
subset <- sample(1:nrow(samples), min(400, nrow(samples)))
layout(matrix(1:2, 1, 2))
matplot(x, potential_curves[, subset], type = "l", lty = 1, 
        col = scales::alpha(1, 1e-2), lwd = 2, ylim = c(-0.2, 0.2),
        main = "potential function", ylab = "", xlab = "population")
lines(x, potential(x = x, a = a, r = r, H = H, Q = Q, K = K), lwd = 2)
matplot(x[-1], dpotential_curves[, subset], type = "l", lty = 1,
        col = scales::alpha(1, 1e-2), lwd = 2, ylim = c(-0.5, 1.3),
        main = "derivative of potential function", ylab = "", xlab = "population")
lines(x[-1], diff(potential(x = x, a = a, r = r, H = H, Q = Q, K = K))/diff(x), lwd = 2)
abline(h = 0, lwd = 2, lty = 3)
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../../figs/scharf_nimble/posterior_deterministic_core_", N_trajectories, "_", seed, ".pdf"), width = 10)
## deterministic core ----
core_curve <- function(a, r, H, Q, K, x0 = 0.3, t = (1:5e3) * 0.5, t.step = 0.5){
  x <- rep(x0, length(t))
  for(t_i in 1:(length(t) - 1)){
    x[t_i + 1] <- x[t_i] + t.step * (x[t_i] * r * (1 - x[t_i] / K)  - a * x[t_i] ^ Q / (x[t_i] ^ Q + H ^ Q))
  }
  return(cbind(t, x))
}
subset <- sample(1:nrow(samples), min(400, nrow(samples)))
deterministic_core <- apply(samples[subset, ], 1, function(row){
  core_curve(a = exp(row['log_a']), r = exp(row['log_r']), 
            H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))[, 2]
})
layout(1)
matplot((1:5e3)/2, deterministic_core, type = "l", lty = 1, 
        col = scales::alpha(1, 1e-1), lwd = 2, 
        ylim = c(0, min(max(deterministic_core, 1.5, na.rm = T), 5)),
        main = "potential function", ylab = "population", xlab = "time")
lines(core_curve(a, r, H, Q, K), lwd = 3, lty = 2, col = "darkgreen")
legend("topleft", lwd = 2, lty = c(1, 2), col = c(scales::alpha(1, 1e-1), "darkgreen"), 
       legend = c("based on posterior draw", "true params"))
## dev.off ----
dev.off()
