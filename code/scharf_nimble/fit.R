## source model.R ----
source("model.R")
## constants ----
constants <- list(
  N = N, N_trajectories = N_trajectories, x0 = x0, 
  t.step = 1/2, N_t = 2*N-1,
  mu_r = log(r), sd_r = 1,
  mu_K = log(K), sd_K = 1,
  mu_a = log(a), sd_a = 1,
  mu_H = log(H), sd_H = 1,
  mu_Q = log(Q), sd_Q = 1,
  mu_sigma = log(sigma), sd_sigma = 1
  # mu_sigma_me = log(sigma_me), sd_sigma_me = 1
)
## inits ----
inits <- list(log_r = log(r), log_K = log(K), 
              log_a = log(a), log_H = log(H), 
              log_Q = log(Q), log_sigma = log(sigma)
              # log_sigma_me = log(sigma_me)
              )
## define model, and compile ----
model <- nimbleModel(code = code, constants = constants, inits = inits)
cmodel <- compileNimble(model)
## set seed + simulate ---- 
seed <- 270
set.seed(seed)
simulate(cmodel, nodes = c(
  # "y", 
  "x", "mu", "sd_x"))
cmodel$setData(
  # "y", 
  "x")
# ## device ----
# pdf(file = paste0("../../figs/scharf_nimble/sim_traj_", N_trajectories, "_", seed, ".pdf"))
# ## plot x ----
# matplot(seq(1, N, l = N/constants$t.step), cmodel$y, 
#         type = "l", ylab = "x", xlab = "time", 
#         col = scales::alpha("black", 0.4), lty = 1)
# matplot(seq(1, N, l = N/constants$t.step), cmodel$x, type = "l",
#         col = scales::alpha("darkred", 0.6), lty = 1, add = T)
# ## dev.off ----
# dev.off()
## specify block sampler ----
mcmcConf <- configureMCMC(cmodel)
head(mcmcConf$getSamplers(), 10)
mcmcConf$removeSamplers()
target <- c("log_r", "log_K", "log_a",
            "log_H", "log_Q", "log_sigma"
            # "log_sigma_me"
)
mcmcConf$addSampler(target = target, type = 'RW_block')
# mcmcConf$addSampler(target = target, type = 'AF_slice')
tail(mcmcConf$getSamplers(), 10)
## compile and run MCMC ----
system.time({
  mcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(mcmc, project = model)
})
n_iterations <- 2e5
system.time({
  Cmcmc$run(n_iterations, nburnin = n_iterations / 2, 
            thin = max(1, n_iterations / 2 / 5e3))
})
samples <- as.matrix(Cmcmc$mvSamples)
## save samples ----
save(samples, file = paste0("data/scharf_nimble/samples_", N_trajectories, "_", seed, ".RData"))
## source plot.R ----
source("code/scharf_nimble/plot.R")