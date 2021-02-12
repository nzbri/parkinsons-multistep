data {
  
  int<lower=0> N;
  real alpha_prior_mean;
  real alpha_prior_sd;
  real C_prior_mean;
  real C_prior_sd;
  real k_prior_mean;
  real k_prior_sd;
  vector[N] x;
  vector[N] y;
  vector[N] y_min;
  vector[N] y_max;
  
  
}
parameters {
  
  real<lower=0> alpha; // exposures
  real<lower=0,upper=100> C; // percent population susceptable
  real<lower=0> k; // steps
  real<lower=0> sigma_sq;
  
}

transformed parameters {
  
  real<lower=0> sigma;
  real<lower=0> y_pred[N];
  
  // Variables for intermediate calcs (for readability, rescaling)
  
  real alpha1;
  real C1;

  sigma = sqrt(sigma_sq);
  
  for (n in 1:N) {
    alpha1 = alpha/1e15;
    C1 = C/100;
    //rates are per 100,000
    y_pred[n] = alpha1 * x[n]^(k-1)*100000/(1+((1-C1)/C1)*exp(alpha1/k*(x[n]^k-1)));
  }
  
  
}

model {
  
  alpha ~ normal(alpha_prior_mean,alpha_prior_sd);
  C ~ normal(C_prior_mean,C_prior_sd);
  k ~ normal(k_prior_mean,k_prior_sd);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  
  y ~ normal(y_pred,sigma);
  
}
