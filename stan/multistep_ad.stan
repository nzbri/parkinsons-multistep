data {
  
  int<lower=0> N;
  real alpha_prior_mean;
  real alpha_prior_sd;
  real k_prior_mean;
  real k_prior_sd;
  vector[N] x;
  vector[N] y;
  vector[N] y_min;
  vector[N] y_max;
  
  
}
parameters {
  real<lower=0,upper=10> alpha; // slope early ages
  real<lower=0> k; // steps
  real<lower=0> sigma_sq;
}

transformed parameters {
  
  real<lower=0> sigma;
  real<lower=0> y_pred[N];
  sigma = sqrt(sigma_sq);
  
  for (n in 1:N) {
    y_pred[n] = (alpha/1e3 * x[n])^(k-1) * 100000;
  }
  
//  print("y_pred = ", y_pred,
//        " alpha = ", alpha,
//        " k = ", k);
  
}

model {
  
  alpha ~ normal(alpha_prior_mean,alpha_prior_sd);
  k ~ normal(k_prior_mean,k_prior_sd);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  
  y ~ normal(y_pred,sigma);
  
}
