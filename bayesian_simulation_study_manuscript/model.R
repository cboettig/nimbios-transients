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
## define stochastic model in BUGS notation ----
code <- nimble::nimbleCode({
  # log(r) ~ dnorm(mu_r, sd_r)
  # log(K) ~ dnorm(mu_K, sd_K)
  # log(a) ~ dnorm(mu_a, sd_a)
  # log(H) ~ dnorm(mu_H, sd_H)
  # log(Q) ~ dnorm(mu_Q, sd_Q)
  # log(sigma) ~ dnorm(mu_sigma, sd_sigma)
  # log(sigma_me) ~ dnorm(mu_sigma_me, sd_sigma_me)
  r ~ dgamma(r_shape, r_rate)
  K ~ dgamma(K_shape, K_rate)
  a ~ dgamma(a_shape, a_rate)
  H ~ dgamma(H_shape, H_rate)
  Q ~ dgamma(Q_shape, Q_rate)
  sigma ~ dgamma(sigma_shape, sigma_rate)
  for(i in 1:N_trajectories){
    x[1, i] <- x0
    for(t in 1:N_t){
      mu[t, i] <- x[t, i] + t.step * (x[t, i] * r * (1 - x[t, i] / K) - 
                                      a * x[t, i] ^ Q / (x[t, i] ^ Q + H ^ Q))
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