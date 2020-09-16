## source simulate.R ----
source("simulate.R")
## source model_functional.R ----
source("model_functional.R")
## source priors and initial values ----
source("constants_inits_functional.R")
## initialize x and create model for fitting (choose INCLUDE_ME) ----
INCLUDE_ME <- TRUE
if(INCLUDE_ME){
inits_functional <- list(
  x = obs_y + rnorm(length(obs_y), sd = diff(range(obs_y)) * 1e-2),
  beta = rep(0, constants_functional$degree)
  # sigma = sigma, sigma_me = sigma_me
)
data <- list(y = obs_y)
} else {
inits_functional <- list(
  beta = rep(0, constants_functional$degree)
  # sigma = sigma, sigma_me = sigma_me
)
data <- list(x = true_x)
}
model_functional <- nimbleModel(code = code_functional, 
                                constants = constants_functional, 
                                inits = inits_functional, data = data)
## compile model ----
cmodel_functional <- compileNimble(model_functional)
## specify samplers and monitors ----
mcmcConf <- configureMCMC(cmodel_functional)
if(INCLUDE_ME) mcmcConf$addMonitors('x')
head(mcmcConf$getSamplers(), 10)
target <- c(paste0("beta[", 1:degree, "]"), "sigma")
if(INCLUDE_ME) target <- c(target, "sigma_me")
mcmcConf$removeSamplers(target)
# mcmcConf$addSampler(target = target, type = 'RW_block')
mcmcConf$addSampler(target = target, type = 'AF_slice')
tail(mcmcConf$getSamplers(), 10)
## compile and run MCMC ----
system.time({
mcmc <- buildMCMC(mcmcConf)
Cmcmc <- compileNimble(mcmc, project = model, resetFunctions = T)
})
n_iterations <- 1e4
system.time({
  Cmcmc$run(n_iterations, nburnin = n_iterations / 2, 
            thin = max(1, n_iterations / 2 / 5e3))
})
samples <- as.matrix(Cmcmc$mvSamples)
## save samples ----
save(samples, file = paste0("data/samples_functional_", N_trajectories, "_", seed, ".RData"))
## source plot.R ----
source("plot.R")