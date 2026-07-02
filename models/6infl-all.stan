data {
  
 // sample size
  int<lower=0> N;//n individuals
  int<lower=0> S;//one or two depending on whether calcification and age effects are sex specific

  // predictors
  array[N] int<lower=1, upper=2> sex;
  vector[N] age_z;
  vector[N] ffm_z;
  vector[N] thoracic_fat_z;
  vector[N] inflam_markers;
}

parameters {
  // intercept
  vector[S] alpha;

  // regression coefficients
  vector[S] gamma;
  real epsilon;
  real eta;

  
  // residual SD
  real<lower=0> sigma;
  
}


model {

  // priors
  alpha ~ normal(0,1);
  gamma ~ normal(0,1);
  epsilon ~ normal(0,1);
  eta ~ normal(0,1);
  sigma ~ exponential(1);
  
  // linear predictor
  vector[N] mu;

  for (i in 1:N) {

    mu[i] =
        alpha[sex[i]]
      + gamma[sex[i]] * age_z[i]
      + epsilon * ffm_z[i]
      + eta * thoracic_fat_z[i];
  }


  // likelihood
  inflam_markers ~ normal(mu, sigma);
}
