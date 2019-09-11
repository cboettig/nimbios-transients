// Fitting mechanistic model (May 1976)
// Fix all parameters but a 

data {
  int<lower=1> N; // how many datapoints
	// vector[N] t; // time, lagged
	vector[N] x; // response
  // vector[N] t1; // time
	}

transformed data{
// p <- list(r = .05, K = 2, Q = 5, H = .38, sigma = .01, a = 0.023, N = 1e4)
   real r = 0.05;
   real K = 2;
   real Q = 5;
   real H = 0.38;
   real sigma = 0.01;
}

parameters {
  real<lower=0> a;   
}

model {	
   real mu[N];
   for(t in 1:(N-1)){
            mu[t] =  x[t] + x[t] * r * (1 - x[t] / K)  - a * x[t] ^ Q / (x[t] ^ Q + H ^ Q);
            x[t+1] ~ lognormal(mu[t], sigma*x[t]); // Sep 2019 -- Carl has sigma*x[t] here (I had just sigma) ... do we need?

   }      
    a ~ normal(0,1); # Can make narrower, but doesn't matter
}

/* generated quantities{
} */
