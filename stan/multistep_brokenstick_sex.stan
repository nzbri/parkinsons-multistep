data {
  
  int<lower=0> N;
  real m_early_prior_mean;
  real m_early_prior_sd;
  real m_late_prior_mean;
  real m_late_prior_sd;
  real bp_prior_mean;
  real bp_prior_sd;
  // If need to be able to modify priors
  //real c_sex_prior_mean; // difference
  //real c_sex_prior_sd; // difference
  //real m_sex_prior_mean; // difference
  //real m_sex_prior_sd; // difference
  vector[N] x;
  vector[N] sex;
  vector[N] y_min;
  vector[N] y_max;
  
  
}
parameters {
  real m_early; // slope early ages
  real m_late; // slope later ages
  real<lower=0.2, upper=0.8> bp; // break point age
  real c;
  real c_sex; // difference in intercept between sexes
  real m_sex; // difference in slope between sexes
  real<lower=0> sigma_sq;
}

transformed parameters {
  
  real<lower=0> sigma;
  real y_pred[N];
  sigma = sqrt(sigma_sq);
  
  for (n in 1:N)
  if (x[n] < bp) {
    y_pred[n] = (m_early + sex[n]*m_sex) * x[n] + 
                c + sex[n] * c_sex;
  } else {
    y_pred[n] = (m_early + sex[n]*m_sex) * bp + 
                (m_late + sex[n]*m_sex) * (x[n]-bp) + 
                c + sex[n] * c_sex;
  }
  
}

model {
  
  m_early ~ normal(m_early_prior_mean,m_early_prior_sd);
  m_late ~ normal(m_late_prior_mean,m_late_prior_sd);
  bp ~ normal(bp_prior_mean,bp_prior_sd);
  c ~ normal(0,40);
  c_sex ~ normal(0,10);
  m_sex ~ normal(0,5);
  sigma_sq ~ normal(0,10);
  
  target += -2 * log(sigma);
  for (n in 1:N)
    target += log(Phi((y_max[n] - y_pred[n]) / sigma)
                  - Phi((y_min[n] - y_pred[n]) / sigma));
}

generated quantities {
  real c_female;
  real m_early_female;
  real m_late_female;
  
  c_female = c + c_sex;
  m_early_female = m_early + m_sex;
  m_late_female = m_late + m_sex;

}
