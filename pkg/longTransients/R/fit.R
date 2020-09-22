#' Create function to fit model using nimble
#'
#' @param data 
#' @param constants 
#' @param inits 
#' @param seed 
#' @param n_iterations in MCMC algorithm
#' @param MVsampler_type either "AF_slice" or "RW_block"
#'
#' @return list including samples from MCMC, simulated data, and timing information
#' @export
fit_mcmc_functional <- function(data, constants, inits, seed, 
                                n_iterations = 1e1, MVsampler_type = "AF_slice"){
  ## initialize x and create model for fitting (INCLUDE_ME) ----
  if(is.null(inits$beta)){
    inits$beta <- rep(0, constants$degree)
  }
  INCLUDE_ME <- !is.null(data$obs_y)
  if(INCLUDE_ME){
    if(is.null(inits$sigma_me) & is.null(inits$x)) inits$sigma_me <- 1e-2
    inits$x <- obs_y + rnorm(length(obs_y), sd = inits$sigma_me)
  }
  code_functional <- make_functional_model_code_nC()
  nimbleModel(code = code_functional, constants = constants, 
              inits = inits, data = data)
  ## compile model ----
  cmodel_functional <- compileNimble(model_functional)
  ## specify samplers and monitors ----
  mcmcConf <- configureMCMC(cmodel_functional)
  if(INCLUDE_ME) mcmcConf$addMonitors('x')
  target <- c(paste0("beta[", 1:degree, "]"), "sigma")
  if(INCLUDE_ME) target <- c(target, "sigma_me")
  mcmcConf$removeSamplers(target)
  mcmcConf$addSampler(target = target, type = MVsampler_type)
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
              sample_time = sample_time, data = data, 
              constants = constants, inits = inits))
}

#' Create function to extract summaries
#'
#' @param samples 
#' @param true 
#' @param xmin 
#' @param xmax 
#'
#' @return scalar; integrated squared error between true derivative of potential and posterior distribution
#' @export
summarize_fit_functional <- function(samples, true, xmin, xmax){
  degree <- length(grep("beta", colnames(samples)))
  x <- seq(xmin, xmax, length.out = 1e2)
  dpotential_curves <- apply(samples[, 1:degree], 1, function(row){
    sapply(x, get_dV, beta = row)
  })
  true_dpotential <- dpotential(x = x, a = true$a, r = true$r, 
                                H = true$H, Q = true$Q, K = true$K)
  ISE <- rowSums((dpotential_curves - true_dpotential)^2) * diff(x[1:2])
}
## ----
fit <- fit_mcmc_functional(seed = 1234, sigma_me = 0.05, N_trajectories = 10, n_iterations = 1e3)
fit$sample_time
ISE <- summarize_fit_functional(samples = fit$samples, true = inits, 
                                xmin = 0.2, xmax = 1.9)
ISE
## save samples ----
save(fit, file = paste0("data/samples_functional_", N_trajectories, "_", seed, ".RData"))
## source plot.R ----
source("plot.R")