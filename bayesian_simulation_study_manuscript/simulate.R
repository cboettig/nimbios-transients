## source model.R ----
source("model.R")
## source priors and initial values ----
source("constants_inits.R")
## define model, and compile ----
INCLUDE_ME <- T
# H <- inits$H <- 0.39
model <- nimbleModel(code = code, constants = constants, inits = inits)
cmodel <- compileNimble(model)
## set seed + simulate ----
seed <- 270
set.seed(seed)
simulate(cmodel, nodes = c("y", "x", "mu", "sd_x"))
obs_y <- cmodel$y
true_x <- cmodel$x
# ## device ----
# pdf(file = paste0("fig/sim_traj_", N_trajectories, "_", seed, ".pdf"))
# plot y, x ----
matplot(seq(1, N, l = N/constants$t.step), obs_y,
        type = "l", ylab = "x", xlab = "time",
        col = scales::alpha("black", 0.4), lty = 1)
matplot(seq(1, N, l = N/constants$t.step), true_x, type = "l",
        col = scales::alpha("darkred", 0.6), lty = 1, add = T)
# ## dev.off ----
# dev.off()
