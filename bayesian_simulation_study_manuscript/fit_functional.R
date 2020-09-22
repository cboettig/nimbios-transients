## create function to fit model using nimble ----
fit_mcmc_functional <- function(seed, sigma_me, N_trajectories, 
                                n_iterations = 1e1, INCLUDE_ME = TRUE){
  ## source simulate.R ----
  source("simulate.R")
  ## source model_functional.R ----
  source("model_functional.R")
  ## source constants_inits_functional.R ----
  source("constants_inits_functional.R")
  ## initialize x and create model for fitting (INCLUDE_ME) ----
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
  compile_time <- system.time({
    mcmc <- buildMCMC(mcmcConf)
    Cmcmc <- compileNimble(mcmc, project = model, resetFunctions = T)
  })
  sample_time <- system.time({
    Cmcmc$run(n_iterations, nburnin = n_iterations / 2, 
              thin = max(1, n_iterations / 2 / 5e3))
  })
  samples <- as.matrix(Cmcmc$mvSamples)
  ## return ----
  return(list(samples = samples, compile_time = compile_time, 
              sample_time = sample_time, obs_y = obs_y))
}
## create function to extract summaries ----
summarize_fit_functional <- function(samples, true, xmin, xmax){
  degree <- length(grep("beta", colnames(samples)))
  x <- seq(xmin, xmax, length.out = 1e2)
  dpotential_curves <- apply(samples[, 1:degree], 1, function(row){
    sapply(x, get_dV, beta = row)
  })
  true_dpotential <- dpotential(x = x, a = true$a, r = true$r, 
                                H = true$H, Q = true$Q, K = true$K)
  ISE <- sum((dpotential_curves - true_dpotential)^2) * diff(x[1:2])
}
## ----
fit <- fit_mcmc_functional(seed = 1234, sigma_me = 0.02, N_trajectories = 10, n_iterations = 1e2)
ISE <- summarize_fit_functional(samples = fit$samples, true = inits, 
                                xmin = 0.2, xmax = 1.9)
ISE
## save samples ----
save(fit, file = paste0("data/samples_functional_", N_trajectories, "_", seed, ".RData"))
## source plot.R ----
source("plot.R")