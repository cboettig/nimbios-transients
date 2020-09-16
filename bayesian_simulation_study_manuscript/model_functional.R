## libraries ----
library(nimble)
## define "spike" normal distribution ----
dspikenorm <- nimbleFunction(
  run = function(x = double(0), mean = double(0), 
                 sd = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if(x < 0){
      if(log) return(-Inf)
      else return(0)
    }
    if(x == 0){
      log_prob <- pnorm(q = 0, mean = mean, sd = sd, log = 1)
      if(log) return(log_prob)
      else return(exp(log_prob))
    }
    log_prob <- dnorm(x = x, mean = mean, sd = sd, log = 1)
    if(log) return(log_prob)
    else return(exp(log_prob))
  })
rspikenorm <- nimbleFunction(
  run = function(n = integer(0, default = 1), mean = double(0),
                 sd = double(0)) {
    returnType(double(0))
    if(n != 1) print("rtruncnorm only allows n = 1; using n = 1.")
    draw <- rnorm(n = 1, mean = mean, sd = sd)
    if(draw <= 0) return(0)
    else return(draw)
  })
get_dV <- nimbleFunction(
  run = function(x = double(0), beta = double(1)) {
    returnType(double(0))
    X <- c(1.524111, 
           7.180e-02 * x - 7.495e-02,
           1.823e-01 * x^2 - 3.805e-01 * x + 1.633e-01,
           4.712e-01 * x^3 - 1.476e+00 * x^2 + 1.376e+00 * x - 3.642e-01,
           1.226 * x^4 - 5.118 * x^3 + 7.402 * x^2 - 4.300 * x + 8.247e-01)
    dV <- inprod(X, beta)
    return(dV)
})
## define stochastic model in BUGS notation ----
code_functional <- nimble::nimbleCode({
  beta[1:degree] ~ dmnorm(mean = beta_mean[1:degree], prec = beta_prec[1:degree, 1:degree])
  sigma ~ dgamma(sigma_shape, sigma_rate)
  for(i in 1:N_trajectories){
    x[1, i] <- x0
    for(t in 1:N_t){
      mu[t, i] <- x[t, i] + t.step * get_dV(x[t, i], beta[1:degree])
      sd_x[t, i] <- sigma * mu[t, i] * sqrt(t.step)
      x[t + 1, i] ~ dspikenorm(mu[t, i], sd_x[t, i])    
    }
  }
  if(INCLUDE_ME){
    sigma_me ~ dgamma(sigma_me_shape, sigma_me_rate)
    for(i in 1:N_trajectories){
      for(t in 1:N){
        y[t, i] ~ dspikenorm(x[t, i], sigma_me)
      }
    }
  }
})