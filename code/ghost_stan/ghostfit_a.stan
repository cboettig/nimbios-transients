// Fitting mechanistic model (May 1976)
// Try to fit all params

data {
  int<lower=1> N; // how many datapoints
  vector[N] x; // response
	}

parameters {
  real<lower=0> r; 
  real<lower=0> K; 
  real<lower=0> Q;   
  real<lower=0> H;   
  real<lower=0> a;   
  real<lower=0> sigma;   
}

/* transformed parameters {

	}*/

model {	
   real mu[N];
   // x[1] = x0;
   for(t in 1:(N-1)){
            mu[t] =  x[t] + x[t] * r * (1 - x[t] / K)  - a * x[t] ^ Q / (x[t] ^ Q + H ^ Q);
            x[t+1] ~ lognormal(mu[t], sigma*x[t]); // Sep 2019 -- Carl has sigma*x[t] here (I had just sigma) ... do we need?
   }
            
    Q ~ normal(5, 0.1);
    H ~ normal(0.38, 0.1);
    a ~ normal(0.0233, 0.1);
    sigma ~ normal(0.01, 0.001);
    // x0 ~ normal(0, 1);

}

/* generated quantities{
} */
