## source simulate.R ----
source("simulate.R")
## source model.R ----
source("model_functional.R")
## source priors and initial values ----
source("constants_inits_functional.R")
## initialize x and create model for fitting (choose INCLUDE_ME) ----
INCLUDE_ME <- TRUE
if(INCLUDE_ME){
  inits <- list(
    # r = r, K = K, a = a, H = H, Q = Q,
    # sigma = sigma, sigma_me = sigma_me
  )
  data <- list(x = true_x)
  model <- nimbleModel(code = code, constants = constants, inits = inits, data = data)
} else {
  inits <- list(
    x = obs_y + rnorm(length(obs_y), sd = diff(range(obs_y)) * 1e-2),
    # r = r, K = K, a = a, H = H, 
    # Q = Q,
    # sigma = sigma
    # sigma_me = sigma_me
  )
  data <- list(y = obs_y)
  model <- nimbleModel(code = code, constants = constants, inits = inits, data = data)
}
## compile model ----
cmodel <- compileNimble(model)
## specify samplers and monitors ----
mcmcConf <- configureMCMC(cmodel)
if(INCLUDE_ME) mcmcConf$addMonitors('x')
# head(mcmcConf$getSamplers(), 10)
# target <- c("r", "K", "a",
#             "H", 
#             "Q",
#             "sigma")
# if(INCLUDE_ME) target <- c(target, "sigma_me")
# mcmcConf$removeSamplers(target)
# # mcmcConf$addSampler(target = target, type = 'RW_block')
# mcmcConf$addSampler(target = target, type = 'AF_slice')
# tail(mcmcConf$getSamplers(), 10)
## compile and run MCMC ----
system.time({
  mcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(mcmc, project = model, resetFunctions = T)
})
  n_iterations <- 1e1
  system.time({
    Cmcmc$run(n_iterations, nburnin = n_iterations / 2, 
              thin = max(1, n_iterations / 2 / 5e3))
  })
  samples <- as.matrix(Cmcmc$mvSamples)
## save samples ----
save(samples, file = paste0("data/samples_functional_", N_trajectories, "_", seed, ".RData"))
## source plot.R ----
source("plot.R")