data {
  int<lower=1> N;
  int<lower=1> S;
  array[N] int<lower=1,upper=S> sex;
  vector[N] age_z;
  vector[N] ffm_z;
  vector[N] thoracic_fat_z;
  array[N] int<lower=0> calc;
}
parameters {
  // intercept
  vector[S] alpha;
  // regression coefficients
  vector[S] gamma;
  real epsilon;
  real eta;
  //zero inflated
  vector[S] alpha_zero;
  vector[S] gamma_zero;
  real epsilon_zero;
  real eta_zero;

  // Lognormal-Poisson mixture: per-observation random effect on log_mu
  // induces extra-Poisson variance (overdispersion)
  real<lower=0> sigma_u;
  vector[N] u_raw; // non-centered parameterization
}
transformed parameters {
  vector[N] u; // realized random effect, u[i] ~ normal(0, sigma_u)
  u = sigma_u * u_raw;
}
model {
  // priors
  alpha ~ normal(0,1);
  gamma ~ normal(0,1);
  epsilon ~ normal(0,1);
  eta ~ normal(0,1);
  //zero inflated
  alpha_zero ~ normal(0,1);
  gamma_zero ~ normal(0,1);
  epsilon_zero ~ normal(0,1);
  eta_zero ~ normal(0,1);


  // overdispersion random effect
  sigma_u ~ normal(0,1); // half-normal due to <lower=0>
  u_raw ~ normal(0,1);   // non-centered: u = sigma_u * u_raw

for(i in 1:N){
  real log_mu;
  real logit_pi;
  log_mu =
      alpha[sex[i]]
    + gamma[sex[i]] * age_z[i]
    + epsilon * ffm_z[i]
    + eta * thoracic_fat_z[i]
    + u[i]; // overdispersion term
  logit_pi =
      alpha_zero[sex[i]]
    + gamma_zero[sex[i]] * age_z[i]
    + epsilon_zero * ffm_z[i]
    + eta_zero * thoracic_fat_z[i];
    
  real pi_i = inv_logit(logit_pi);
  
  if(calc[i] == 0){
    target += log_sum_exp(
      bernoulli_lpmf(1 | pi_i),
      bernoulli_lpmf(0 | pi_i)
      + poisson_log_lpmf(0 | log_mu)
    );
  } else {
    target +=
      bernoulli_lpmf(0 | pi_i)
      + poisson_log_lpmf(calc[i] | log_mu);
  }
}
}
