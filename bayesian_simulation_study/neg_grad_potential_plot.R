## libraries ----
library(longTransients)
## posteriors ------------------------------------------------------------------
as <- c(0.0225, 0.023, 0.0235)
subsets <- list(1, 1:2, 1:5, 1:10)
n_trajs <- c(1, 2, 5, 10)
samples <- lapply(as, function(a){
  lapply(subsets, function(subset){
    load(paste0("data/fit_parametric_sigma_me_0_01_a_", gsub(".", "_", a, fixed = T), 
                "_y_subset_", paste(subset, collapse = "_"), ".RData"))
    return(fit)
  })
})
## compare_potential curves ----
## device ----
pdf("fig/potential_curves_a.pdf", height = 6, width = 7)
## plot ----
x_eval <- seq(0.2, 1.8, l = 2e2)
true <- list(r = 0.05, K = 2, H = 0.38, Q = 5, sigma_me = 0.01)
layout(matrix(1:9, 3, 3))
par(oma = c(3, 6, 6, 0), mar = c(1, 1, 0, 0.1))
for(n_traj in n_trajs[-3]){
  for(a in as){
    sample_a_subset <- samples[[which(a == as)]][[which(n_traj == n_trajs)]]
    true$a <- a
    xlab <- ""; xaxt <- "n"; yaxt <- "n"
    if(a == as[3]) {xlab <- "population"; xaxt <- "s"}
    if(n_traj == 1) yaxt <- "s"
    plot_potential(samples = sample_a_subset$samples, true = true, x = x_eval, xlab = xlab,
                   ylim_grad = 0.02 * c(-1, 1), which_plot = "grad", main = "", xaxt = xaxt, yaxt = yaxt)
    if(n_traj == 1) mtext(a, 2, line = 3)
    if(a == as[1]) mtext(n_traj, 3, line = 1)
    if(a == as[3]) mtext("population", 1, line = 2.8)
  }
}
title(main = "number of trajectories", outer = T, cex.main = 2, line = 4)
title(ylab = "a", outer = T, line = 4, xpd = T, cex.lab = 2.5)
## dev.off ----
dev.off()