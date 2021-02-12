data {
  
  int<lower=0> N;
  real alpha_prior_mean;
  real alpha_prior_sd;
  real beta_prior_mean;
  real beta_prior_sd;
  real k_prior_mean;
  real k_prior_sd;
  vector[N] x;
  vector[N] y;
  vector[N] y_min;
  vector[N] y_max;
  
  
}
parameters {
  real<lower=0,upper=0.01> alpha; // slope early ages
  real<lower=0,upper=0.02> beta; // drop off at older ages
  real<lower=0> k; // steps
  real<lower=0> sigma_sq;
}

transformed parameters {
  
  real<lower=0> sigma;
  // y_pred[N] should ideally be <lower=0> but for the beta model it is mathematically
  // possible to cross zero at older ages and don't want these rejected
  real y_pred[N]; 

  
  sigma = sqrt(sigma_sq);
  
  for (n in 1:N) {
    y_pred[n] = (alpha * x[n])^(k-1) * (1-beta*x[n]) * 100000;
  }
  
  //print("y_pred = ", y_pred,
  //      " alpha = ", alpha,
  //      " beta = ", beta,
  //      " k = ", k);
  
}

model {
  
  alpha ~ normal(alpha_prior_mean,alpha_prior_sd);
  beta ~ normal(beta_prior_mean,beta_prior_sd);
  k ~ normal(k_prior_mean,k_prior_sd);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  
  y ~ normal(y_pred,sigma);
  
}
