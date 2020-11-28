#' Fit functional model using nimble
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
fit_mcmc_functional <- function(data, constants, inits = NULL, seed, 
                                n_iterations = 1e1, MVsampler_type = "AF_slice"){
  ## initialize x and create model for fitting (INCLUDE_ME) ----
  if(is.null(inits$beta)){
    inits$beta <- rep(0, constants$degree)
  }
  INCLUDE_ME <- !is.null(data$y)
  if(INCLUDE_ME){
    if(is.null(inits$sigma_me) & is.null(inits$x)) inits$sigma_me <- 1e-2
    inits$x <- data$y + rnorm(length(data$y), sd = inits$sigma_me)
  }
  code_functional <- make_model_code_functional_nC()
  model_functional <- nimbleModel(code = code_functional, constants = constants, 
                                  inits = inits, data = data)
  ## compile model ----
  cmodel_functional <- compileNimble(model_functional)
  ## specify samplers and monitors ----
  mcmcConf <- configureMCMC(cmodel_functional)
  if(INCLUDE_ME) mcmcConf$addMonitors('x')
  target <- c(paste0("beta[", 1:constants$degree, "]"), "sigma")
  if(INCLUDE_ME) target <- c(target, "sigma_me")
  mcmcConf$removeSamplers(target)
  mcmcConf$addSampler(target = target, type = MVsampler_type)
  ## compile and run MCMC ----
  compile_time <- system.time({
    mcmc <- buildMCMC(mcmcConf)
    Cmcmc <- compileNimble(mcmc, project = model_functional, resetFunctions = T)
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

#' Fit parametric model using nimble
#'
#' @param data 
#' @param constants 
#' @param inits 
#' @param seed 
#' @param n_iterations 
#' @param MVsampler_type 
#'
#' @return
#' @export
fit_mcmc <- function(data, constants, inits = NULL, seed, 
                     n_iterations = 1e1, MVsampler_type = "AF_slice"){
  ## initialize x and create model for fitting (INCLUDE_ME) ----
  INCLUDE_ME <- !is.null(data$y)
  if(INCLUDE_ME){
    if(is.null(inits$sigma_me) & is.null(inits$x)) inits$sigma_me <- 1e-2
    inits$x <- data$y + rnorm(length(data$y), sd = inits$sigma_me)
  }
  code <- make_model_code_nC()
  model <- nimbleModel(code = code, constants = constants, 
                       inits = inits, data = data)
  ## compile model ----
  cmodel <- compileNimble(model)
  ## specify samplers and monitors ----
  mcmcConf <- configureMCMC(cmodel)
  if(INCLUDE_ME) mcmcConf$addMonitors('x')
  target <- c("r", "K", "a", "H", "Q", "sigma")
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

#' Growth function
#'
#' @param x 
#' @param r 
#' @param K 
#'
#' @return
growth <- function(x, r, K){x * r * (1 - x / K)}

#' Consumption function
#'
#' @param x 
#' @param a 
#' @param H 
#' @param Q 
#'
#' @return
consumption <- function(x, a, H, Q){a * x^Q / (x^Q + H^Q)}

#' Compute derivative of potential function 
#'
#' @param x 
#' @param a 
#' @param r 
#' @param H 
#' @param Q 
#' @param K 
#'
#' @return
#' @export
dpotential <- function(x = seq(0, 2, l = 1e2), a, r, H, Q, K){
  growth(x = x, r = r, K = K) - 
    consumption(x = x, a = a, H = H, Q = Q)
}

#' Compute potential function by numerically integrating derivative calculated with dpotential()
#'
#' @param x 
#' @param a 
#' @param r 
#' @param H 
#' @param Q 
#' @param K 
#'
#' @return
#' @export
potential <- function(x = seq(0, 2, l = 1e2), a, r, H, Q, K){
  out <- cumsum(-dpotential(x = x, r = r, K = K, a = a, H = H, Q = Q)[-1] * diff(x))
  out <- out - out[1]
}

#' Extract ISE for parametric model
#'
#' @param samples 
#' @param true 
#' @param x 
#'
#' @return scalar; integrated squared error between true derivative of potential and posterior distribution
#' @export
get_ISE <- function(samples, true, x = seq(0, 2, length.out = 1e2)){
  dpotential_curves <- apply(samples, 1, function(row){
    sapply(x, dpotential, a = row['a'], r = row['r'], 
           H = row['H'], Q = row['Q'], K = row['K'])
  })
  true_dpotential <- dpotential(x = x, a = true$a, r = true$r, 
                                H = true$H, Q = true$Q, K = true$K)
  ISE <- rowSums((dpotential_curves - true_dpotential)^2) * diff(x[1:2])
}

#' Extract ISE for functional model
#'
#' @param samples 
#' @param true 
#' @param x
#'
#' @return
#' @export
get_ISE_functional <- function(samples, true, x = seq(0, 2, length.out = 1e2)){
  degree <- length(grep("beta", colnames(samples)))
  dpotential_curves <- apply(samples[, 1:degree], 1, function(row){
    sapply(x, dV, beta = row)
  })
  true_dpotential <- dpotential(x = x, a = true$a, r = true$r, 
                                H = true$H, Q = true$Q, K = true$K)
  ISE <- rowSums((dpotential_curves - true_dpotential)^2) * diff(x[1:2])
}

#' Compute summaries of effective sample size for derivative of potential for parametric model
#'
#' @param samples 
#' @param x 
#' @param na.rm default `TRUE`
#' @importFrom coda effectiveSize
#'
#' @return vector of min, mean, and median ESS across derivative of potential
#' @export
get_ESS <- function(samples, x = seq(0.1, 1.9, length.out = 1e2), na.rm = TRUE){
  dpotential_curves <- apply(samples, 1, function(row){
    sapply(x, dpotential, a = row['a'], r = row['r'], 
           H = row['H'], Q = row['Q'], K = row['K'])
  })
  ESS <- effectiveSize(t(dpotential_curves))
  return(c(min = min(ESS, na.rm = na.rm), 
           mean = mean(ESS, na.rm = na.rm), 
           median = median(ESS, na.rm = na.rm)))
}

#' Compute summaries of effective sample size for derivative of potential
#'
#' @param samples 
#' @param x 
#'
#' @return vector of min, mean, and median ESS across derivative of potential for functional model
#' @export
get_ESS_functional <- function(samples, x = seq(0.1, 1.9, length.out = 1e2)){
  degree <- length(grep("beta", colnames(samples)))
  dpotential_curves <- apply(samples[, 1:degree], 1, function(row){
    sapply(x, dV, beta = row)
  })
  ESS <- effectiveSize(t(dpotential_curves))
  return(c(min = min(ESS), mean = mean(ESS), median = median(ESS)))
}
