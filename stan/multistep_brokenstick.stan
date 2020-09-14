data {
  
  int<lower=0> N;
  real m_early_prior_mean;
  real m_early_prior_sd;
  real m_late_prior_mean;
  real m_late_prior_sd;
  real bp_prior_mean;
  real bp_prior_sd;
  vector[N] x;
  vector[N] y_min;
  vector[N] y_max;
  
  
}
parameters {
  real m_early; // slope early ages
  real m_late; // slope later ages
  real<lower=0.2, upper=0.8> bp; // break point age
  real c;
  real<lower=0> sigma_sq;
}

transformed parameters {
  
  real<lower=0> sigma;
  real y_pred[N];
  sigma = sqrt(sigma_sq);
  
  for (n in 1:N)
  if (x[n] < bp) {
    y_pred[n] = m_early * x[n] + c;
  } else {
    y_pred[n] = m_early * bp + m_late * (x[n]-bp) + c;
  }
  
}

model {
  
  m_early ~ normal(m_early_prior_mean,m_early_prior_sd);
  m_late ~ normal(m_late_prior_mean,m_late_prior_sd);
  bp ~ normal(bp_prior_mean,bp_prior_sd);
  c ~ normal(0,40);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  for (n in 1:N)
    target += log(Phi((y_max[n] - y_pred[n]) / sigma)
                  - Phi((y_min[n] - y_pred[n]) / sigma));
}
