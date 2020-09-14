data {
  int<lower=0> N;
  real m_prior_mean;
  real m_prior_sd;
  vector[N] x;
  vector[N] y_min;
  vector[N] y_max;
  
  
}
parameters {
  real m;
  real c;
  real<lower=0> sigma_sq;
}

transformed parameters {
  real<lower=0> sigma;
  real y_pred[N];
  sigma = sqrt(sigma_sq);
  
  for (n in 1:N)
    y_pred[n] = m * x[n] + c;
  
}

model {
  
  m ~ normal(m_prior_mean,m_prior_sd);
  c ~ normal(0,40);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  for (n in 1:N)
    target += log(Phi((y_max[n] - y_pred[n]) / sigma)
                  - Phi((y_min[n] - y_pred[n]) / sigma));
}
