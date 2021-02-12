data {
  
  int<lower=0> N;
  int<lower=0> DEBUG;
  real alpha_prior_mean;
  real alpha_prior_sd;
  real C_prior_mean;
  real C_prior_sd;
  real k_prior_mean;
  real k_prior_sd;
  vector[N] x;
  vector[N] y;
  vector[N] sex;
  vector[N] y_min;
  vector[N] y_max;
  
}

parameters {
  
  real<lower=0> alpha; // product of exposures
  real alpha_sex; // product of exposures (sex effect)
  real<lower=0,upper=100> C; // percent population susceptable
  real C_sex; // population susceptable (sex effect)
  real<lower=0> k; // steps
  real k_sex; // steps (sex efect)
  real<lower=0> sigma_sq;
  
}

transformed parameters {
  
  real<lower=0> sigma;
  real<lower=0> y_pred[N];
  
  // Variables for intermediate calcs (for readability, rescaling)
  
  real alpha1;
  real k1;
  real C1;
  real numerator;
  real denominator; 
  
  sigma = sqrt(sigma_sq);
  
  for (n in 1:N) {
    
    alpha1 = (alpha + alpha_sex*sex[n])/1e15;
    k1 = k+k_sex*sex[n];
    C1 = (C+C_sex*sex[n])/100;
    
    numerator = alpha1 * x[n]^(k1-1)*100000; // incidence is per 100000
    denominator = (1+((1-C1)/C1) * exp(alpha1/k1*(x[n]^k1-1)));
    
    y_pred[n] = numerator / denominator;
                
  }
  
  if (DEBUG) {
    
    print("y_pred = ", y_pred,
          " alpha = ", alpha,
          " C = ", C,
          " k = ", k);
  }
  
}

model {
  
  alpha ~ normal(alpha_prior_mean,alpha_prior_sd);
  C ~ normal(C_prior_mean,C_prior_sd);
  k ~ normal(k_prior_mean,k_prior_sd);
  alpha_sex ~ normal(0,10);
  C_sex ~ normal(0,5);
  k_sex ~ normal(0,1);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  
  y ~ normal(y_pred,sigma);
  
}

generated quantities {
  real alpha_female;
  real C_female;
  real k_female;
  
  alpha_female = alpha + alpha_sex;
  C_female = C + C_sex;
  k_female = k + k_sex;

}
