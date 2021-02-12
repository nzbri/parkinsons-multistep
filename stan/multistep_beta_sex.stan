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
  vector[N] sex;
  
}
parameters {
  real<lower=0,upper=10> alpha; // slope early ages
  real<lower=0> k; // steps
  real<lower=0,upper=3> beta; // drop off at older ages
  real<lower=0> sigma_sq;
  real alpha_sex; // product of exposures (sex effect)
  real k_sex; // steps (sex efect)
  real<lower=0, upper=3> beta_sex; // drop off at older ages (sex effect)
}

transformed parameters {
  
  real<lower=0> sigma;
  real<lower=0> y_pred[N];
  
  // Variables for intermediate calcs (for readability, rescaling)
  
  real alpha1;
  real k1;
  real beta1;
  
  sigma = sqrt(sigma_sq);  
  
  for (n in 1:N) {
    
    alpha1 = (alpha + alpha_sex*sex[n])/1e3;
    k1 = k+k_sex*sex[n];
    beta1 = (beta+beta_sex*sex[n])/1e2;
    
    y_pred[n] = (alpha1 * x[n])^(k1-1) * (1-beta1*x[n]) * 100000;
  }
  
}

model {
  
  alpha ~ normal(alpha_prior_mean,alpha_prior_sd);
  k ~ normal(k_prior_mean,k_prior_sd);
  beta ~ normal(beta_prior_mean,beta_prior_sd);
  alpha_sex ~ normal(0,3);
  k_sex ~ normal(0,3);
  beta_sex ~ normal(0,1);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  
  y ~ normal(y_pred,sigma);
  
}



generated quantities {
  
  real alpha_female;
  real k_female;
  real beta_female;
  
  alpha_female = alpha + alpha_sex;
  k_female = k + k_sex;
  beta_female = beta + beta_sex;

}
