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

  real<lower=0> phi;
  
}

model {

  vector[N] log_mu;

  // priors
  alpha ~ normal(0,1);
  gamma ~ normal(0,1);
  epsilon ~ normal(0,1);
  eta ~ normal(0,1);

  phi ~ exponential(1);

  for(i in 1:N){

    log_mu[i] =

        alpha[sex[i]]
      + gamma[sex[i]] * age_z[i]
      + epsilon * ffm_z[i]
      + eta * thoracic_fat_z[i];

  }

  calc ~ neg_binomial_2_log(log_mu, phi);
}

