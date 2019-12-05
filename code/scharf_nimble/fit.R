## source model.R ----
source("model.R")
## constants ----
constants <- list(
  N = N, x0 = x0,
  mu_r = log(r), sd_r = 1,
  mu_K = log(K), sd_K = 1,
  mu_a = log(a), sd_a = 1,
  mu_H = log(H), sd_H = 1,
  mu_Q = log(Q), sd_Q = 1,
  mu_sigma = log(sigma), sd_sigma = 1
)
## inits ----
inits <- list(log_r = log(r), log_K = log(K), 
              log_a = log(a), log_H = log(H), 
              log_Q = log(Q), log_sigma = log(sigma))
## define model, set data, and compile ----
model <- nimbleModel(code = code, constants = constants, inits = inits)
cmodel <- compileNimble(model)
## set seed + simulate ----
seed <- 270
set.seed(seed)
simulate(cmodel, nodes = c("x", "mu", "sd_x"))
cmodel$setData("x")
## specify block sampler ----
mcmcConf <- configureMCMC(cmodel)
mcmcConf$getSamplers()
mcmcConf$removeSamplers(c("log_r", "log_K", "log_a",
                          "log_H", "log_Q", "log_sigma"))
mcmcConf$addSampler(target = c("log_r", "log_K", "log_a",
                               "log_H", "log_Q", "log_sigma"),
                    type = 'RW_block')
mcmcConf$getSamplers()
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
save(samples, file = paste0("../../data/scharf_nimble/samples_", seed, ".RData"))
## source plot.R ----
source("plot.R")