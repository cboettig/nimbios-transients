## compare n = 1 and n = 10 ----
## load samples ----
source("model.R")
seed <- 270
load(paste0("../../data/scharf_nimble/samples_", 1, "_", seed, ".RData"))
samples_1 <- samples
load(paste0("../../data/scharf_nimble/samples_", 10, "_", seed, ".RData"))
samples_10 <- samples
rm(samples); gc()
## device ----
pdf(file = paste0("../../figs/scharf_nimble/compare_plots_", seed, ".pdf"),
    height = 6)
## densities ----
colors <- RColorBrewer::brewer.pal(5, "Greys")[3:5]
ltys <- 1:3
names(colors) <- names(ltys) <- c('prior', '1', '10')
## layout 
layout(matrix(1:6, 3, 2))
par(mar = c(3, 3, 2, 1))
## r
r_seq <- seq(0, 0.17, l = 1e2)
ylims <- range(density(exp(samples_1[, 'log_r']))$y,
               density(exp(samples_10[, 'log_r']))$y)
plot(r_seq, dlnorm(r_seq, log(r), 1), type = "l", col = colors['prior'], 
     lwd = 2, lty = ltys['prior'], ylim = ylims,
     main = "r", ylab = "density")
lines(density(exp(samples_1[, 'log_r'])), col = colors['1'], 
      lwd = 2, lty = ltys['1'])
lines(density(exp(samples_10[, 'log_r'])), col = colors['10'], 
      lwd = 2, lty = ltys['10'])
legend("topright", lty = ltys, col = colors, lwd = 2, 
       legend = names(colors), bty = "n")
points(r, par()$usr[3], lwd = 2, xpd = T, pch = 4)
## K
K_seq <- seq(0, 3.5, l = 1e2)
ylims <- range(density(exp(samples_1[, 'log_K']))$y,
               density(exp(samples_10[, 'log_K']))$y)
plot(K_seq, dlnorm(K_seq, log(K), 1), type = "l", col = colors['prior'], 
     lwd = 2, lty = ltys['prior'], ylim = ylims,
     main = "K", ylab = "density")
lines(density(exp(samples_1[, 'log_K'])), col = colors['1'], 
      lwd = 2, lty = ltys['1'])
lines(density(exp(samples_10[, 'log_K'])), col = colors['10'], 
      lwd = 2, lty = ltys['10'])
points(K, par()$usr[3], lwd = 2, xpd = T, pch = 4)
## a
a_seq <- seq(0, 0.08, l = 1e2)
ylims <- range(density(exp(samples_1[, 'log_a']))$y,
               density(exp(samples_10[, 'log_a']))$y)
plot(a_seq, dlnorm(a_seq, log(a), 1), type = "l", col = colors['prior'], 
     lwd = 2, lty = ltys['prior'], ylim = ylims,
     main = "a", ylab = "density")
lines(density(exp(samples_1[, 'log_a'])), col = colors['1'], 
      lwd = 2, lty = ltys['1'])
lines(density(exp(samples_10[, 'log_a'])), col = colors['10'], 
      lwd = 2, lty = ltys['10'])
points(a, par()$usr[3], lwd = 2, xpd = T, pch = 4)
## H
H_seq <- seq(0.25, 0.5, l = 1e2)
ylims <- range(density(exp(samples_1[, 'log_H']))$y,
               density(exp(samples_10[, 'log_H']))$y)
plot(H_seq, dlnorm(H_seq, log(H), 1), type = "l", col = colors['prior'], 
     lwd = 2, lty = ltys['prior'], ylim = ylims,
     main = "h", ylab = "density")
lines(density(exp(samples_1[, 'log_H'])), col = colors['1'], 
      lwd = 2, lty = ltys['1'])
lines(density(exp(samples_10[, 'log_H'])), col = colors['10'], 
      lwd = 2, lty = ltys['10'])
points(H, par()$usr[3], lwd = 2, xpd = T, pch = 4)
## Q
Q_seq <- seq(0, 15, l = 1e2)
ylims <- range(density(exp(samples_1[, 'log_Q']))$y,
               density(exp(samples_10[, 'log_Q']))$y)
plot(Q_seq, dlnorm(Q_seq, log(Q), 1), type = "l", col = colors['prior'], 
     lwd = 2, lty = ltys['prior'], ylim = ylims,
     main = "q", ylab = "density")
lines(density(exp(samples_1[, 'log_Q'])), col = colors['1'], 
      lwd = 2, lty = ltys['1'])
lines(density(exp(samples_10[, 'log_Q'])), col = colors['10'], 
      lwd = 2, lty = ltys['10'])
points(Q, par()$usr[3], lwd = 2, xpd = T, pch = 4)
## sigma
sigma_seq <- seq(0.02 - 1e-3, 0.02 + 1e-3, l = 1e2)
ylims <- range(density(exp(samples_1[, 'log_sigma']))$y,
               density(exp(samples_10[, 'log_sigma']))$y)
plot(sigma_seq, dlnorm(sigma_seq, log(sigma), 1), type = "l", col = colors['prior'], 
     lwd = 2, lty = ltys['prior'], ylim = ylims,
     main = expression(sigma), ylab = "density")
lines(density(exp(samples_1[, 'log_sigma'])), col = colors['1'], 
      lwd = 2, lty = ltys['1'])
lines(density(exp(samples_10[, 'log_sigma'])), col = colors['10'], 
      lwd = 2, lty = ltys['10'])
points(sigma, par()$usr[3], lwd = 2, xpd = T, pch = 4)

## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../../figs/scharf_nimble/compare_posterior_potentials_", seed, ".pdf"))
## calc potential curves ----
growth <- function(x, r, K){x * r * (1 - x / K)}
consumption <- function(x, a, H, Q){a * x^Q / (x^Q + H^Q)}
x <- seq(0, 2, l=1e2)
potential <- function(x = seq(0, 2, l = 1e2), a, r, H, Q, K){
  - cumsum(growth(x = x, r = r, K = K) - 
             consumption(x = x, a = a, H = H, Q = Q))
}
potential_curves_1 <- apply(samples_1, 1, function(row){
  potential(x = x, a = exp(row['log_a']), r = exp(row['log_r']), 
            H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))
})
dpotential_curves_1 <- apply(samples_1, 1, function(row){
  p <- potential(x = x, a = exp(row['log_a']), r = exp(row['log_r']), 
                 H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))
  diff(p) / diff(x)
})
potential_curves_10 <- apply(samples_10, 1, function(row){
  potential(x = x, a = exp(row['log_a']), r = exp(row['log_r']), 
            H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))
})
dpotential_curves_10 <- apply(samples_10, 1, function(row){
  p <- potential(x = x, a = exp(row['log_a']), r = exp(row['log_r']), 
                 H = exp(row['log_H']), Q = exp(row['log_Q']), K = exp(row['log_K']))
  diff(p) / diff(x)
})
## plot potential curves ----
subset <- sample(1:nrow(samples_1), min(400, nrow(samples_1)))
layout(matrix(1:4, 2, 2), width = c(1, 0.82))
par(mar = c(1, 4, 3, 1))
matplot(x, potential_curves_1[, subset], type = "l", lty = 1, 
        col = scales::alpha("black", 1e-2), lwd = 2, ylim = c(-0.2, 0.2),
        ylab = "potential", main = "1 realization", xlab = "", xaxt = "n")
lines(x, potential(x = x, a = a, r = r, H = H, Q = Q, K = K), lwd = 2)
par(mar = c(4, 4, 0, 1))
matplot(x[-1], dpotential_curves_1[, subset], type = "l", lty = 1,
        col = scales::alpha("black", 1e-2), lwd = 2, ylim = c(-0.5, 1.3),
        ylab = "derivative of potential", xlab = "population")
lines(x[-1], diff(potential(x = x, a = a, r = r, H = H, Q = Q, K = K))/diff(x), lwd = 2)
abline(h = 0, lwd = 2, lty = 3)
par(mar = c(1, 0, 3, 1))
matplot(x, potential_curves_10[, subset], type = "l", lty = 1, 
        col = scales::alpha("black", 1e-2), lwd = 2, ylim = c(-0.2, 0.2),
        ylab = "", main = "10 realizations", xlab = "", xaxt = "n", yaxt = "n")
lines(x, potential(x = x, a = a, r = r, H = H, Q = Q, K = K), lwd = 2)
par(mar = c(4, 0, 0, 1))
matplot(x[-1], dpotential_curves_10[, subset], type = "l", lty = 1,
        col = scales::alpha("black", 1e-2), lwd = 2, ylim = c(-0.5, 1.3),
        ylab = "", xlab = "population", yaxt = "n")
lines(x[-1], diff(potential(x = x, a = a, r = r, H = H, Q = Q, K = K))/diff(x), lwd = 2)
abline(h = 0, lwd = 2, lty = 3)
## dev.off ----
dev.off()
