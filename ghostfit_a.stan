// Fitting mechanistic model (May 1976)
// Note, this is not good ... 
// To do: I set r and K impressively hack-y ... and fix x0
// Tonight, fix all values (in param block) but a ... 

data {
  int<lower=1> N; // how many datapoints
	// vector[N] t; // time, lagged
	vector[N] x; // response
  // vector[N] t1; // time
	}

parameters {
  real<lower=0.049, upper=0.051> r; // r = 0.05 # lognormal(0.05, .01)
  real<lower=1.9, upper=2.1> K; //  K = 2 # lognormal(2, .01)
  real<lower=4.9, upper=5.1> Q;   
  real<lower=0.37, upper=0.39> H;   
  real<lower=0> a;   
  real<lower=0> x0;
  real<lower=0> sigma;   
}

/* transformed parameters {

	}*/

model {	
   real mu[N];
   // x[1] = x0;
   for(t in 1:(N-1)){
            mu[t] =  x[t] + x[t] * r * (1 - x[t] / K)  - a * x[t] ^ Q / (x[t] ^ Q + H ^ Q);
            x[t+1] ~ lognormal(mu[t], sigma);
   }
            
    Q ~ normal(5, 0.1);
    H ~ normal(0.38, 0.1);
    a ~ normal(0.0233, 0.1);
    sigma ~ normal(0.01, 0.001);
    x0 ~ normal(0, 1);

}

/* generated quantities{
} */
