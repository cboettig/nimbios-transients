## libraries ----
library(FNN)
library(xtable)
## posteriors ------------------------------------------------------------------
load("data/fit_parametric_sigma_me_0_01_a_0_023_y_subset_1.RData")
samples_1 <- fit$samples
load("data/fit_parametric_sigma_me_0_01_a_0_023_y_subset_1_2.RData")
samples_2 <- fit$samples
load("data/fit_parametric_sigma_me_0_01_a_0_023_y_subset_1_2_3_4_5.RData")
samples_5 <- fit$samples
load("data/fit_parametric_sigma_me_0_01_a_0_023_y_subset_1_2_3_4_5_6_7_8_9_10.RData")
samples_10 <- fit$samples
rm(fit)
nmcmc <- nrow(samples_1)

## KL --------------------------------------------------------------------------
k <- 10 # KL divergence max number of nearest neighbors to search
par_nm <- c("r", "K", "a", "H", "Q", "sigma", "sigma_me")
n_trajs <- c(1, 2, 5, 10)
kl_tab <- sapply(n_trajs, function(n_traj){
  sapply(par_nm, function(par){
    prior <- rgamma(nmcmc, as.numeric(constants_fit[paste0(par, "_shape")]), 
                    as.numeric(constants_fit[paste0(par, "_rate")]))
    KL.divergence(prior, get(paste0("samples_", n_traj))[, par], k)[k]
  })
})

kl_tab

xtable(kl_tab) # create latex table

## ----
## marginal posterior densities  ----------------------------------------------------------------------
## device ----
pdf("fig/KL_divergence.pdf", height = 9)
## plot ----
n_traj_colors <- rev(gray.colors(5))[-1]
prior_col <- gray.colors(5)[5]
layout(matrix(1:8, 4, 2))
par(mar = c(3, 4, 3, 1))
par_ranges <- data.frame("r" = c(0, 0.28), "a" = c(0, 0.27), "Q" = c(0, 17), 
                         "K" = c(1.2, 4), "H" = c(0.33, 0.55), "sigma" = c(0.018, 0.0225),
                         "sigma_me" = c(0.007, 0.012))
for(param in par_nm){
  post_densities <- sapply(n_trajs, function(n_traj){
    par_sample <- get(paste0("samples_", n_traj))[, param]
    density(par_sample)
  }, simplify = F)
  post_densities_range <- range(lapply(post_densities, function(d) range(d$y)))
  par_range <- par_ranges[, param]
  par_seq <- seq(par_range[1], par_range[2], l = 2e2)
  prior_density <- dgamma(par_seq, as.numeric(constants_fit[paste0(param, "_shape")]), 
                          as.numeric(constants_fit[paste0(param, "_rate")]))
  main <- param
  if(param == "Q") main <- "q"
  if(param == "H") main <- "h"
  if(length(grep("sigma", param)) > 0) main <- expression(sigma)
  if(length(grep("sigma_me", param)) > 0) main <- expression(sigma[me])
  plot(par_seq, prior_density, type = "l", lty = 1, lwd = 2, col = prior_col,
       ylab = "Density", xlab = "", ylim = post_densities_range, main = main)
  sapply(1:length(n_trajs), function(i){
    lines(post_densities[[i]], col = n_traj_colors[i], lwd = 2, lty = i + 1)
  })
  points(inits_sim[param], par()$usr[3], xpd = T, pch = 4, lwd = 2)
  # if(param == "r"){
  #   legend("top", col = n_traj_colors, lwd = 2, lty = 2:5,
  #          legend = n_trajs, title = "no. trajectories", bty = "n")
  # }
  legend("topright", legend = c("prior", paste0(round(kl_tab[param, ], 3), " (", n_trajs, ")")), 
         bty = "n", title = "KL div. (no. traj.)", lwd = 2, lty = 1:5, col = c(prior_col, n_traj_colors))
}
## dev.off ----
dev.off()
## ----
## compare_potential curves ----
## device ----
pdf("fig/potential_curves.pdf", height = 6, width = 10)
## plot ----
x_eval <- seq(0.2, 1.8, l = 2e2)
layout(matrix(1:8, 2, 4))
par(oma = c(0, 0.7, 0, 0), mar = c(3, 4, 1, 1))
for(n_traj in n_trajs){
  plot_potential(samples = get(paste0("samples_", n_traj)), true = inits_sim, x = x_eval,
                 ylim_potential = c(-0.005, 0.002), which_plot = "potential", 
                 main = paste(n_traj, "realizations"))
  if(n_traj == 1) mtext("potential", 2, 3)
  plot_potential(samples = get(paste0("samples_", n_traj)), true = inits_sim, x = x_eval,
                 ylim_grad = 0.02 * c(-1, 1), which_plot = "grad", main = "")
  if(n_traj == 1) mtext("derivative of potential", 2, 3)
}
## dev.off ----
dev.off()