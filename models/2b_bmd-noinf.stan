data {
  
 // sample size
  int<lower=0> N;//n individuals
  int<lower=0> C;//n calcification markers
  int<lower=0> S;//one or two depending on whether calcification and age effects are sex specific

  // outcome
  vector[N] bmd_z;

  // predictors
  array[N] int<lower=1, upper=2> sex;
  vector[N] age_z;
  vector[N] ffm_z;
  vector[N] thoracic_fat_z;
  
  // calcification indicators
  matrix[N, C] calc_markers; // 1= cac, 2 = tac, 3 = aac
  vector[N] any_calc;
  
  }

parameters {
  // intercept
  vector[S] alpha;

  // regression coefficients
  vector[S] beta1;
  vector[S] beta2;
  vector[S] gamma;
  real epsilon;
  real eta;

  
  // residual SD
  real<lower=0> sigma;
  
}

transformed parameters {

  vector[N] calc;
  calc = calc_markers * rep_vector(1.0 / C, C);

}

model {

  // priors
  alpha ~ normal(0,1);
  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  gamma ~ normal(0,1);
  epsilon ~ normal(0,1);
  eta ~ normal(0,1);
  sigma ~ exponential(1);
  
  // linear predictor
  vector[N] mu;

  for (i in 1:N) {

    mu[i] =
        alpha[sex[i]]
      + beta1[sex[i]] * any_calc[i]
      + beta2[sex[i]] * calc[i]
      + gamma[sex[i]] * age_z[i]
      + epsilon * ffm_z[i]
      + eta * thoracic_fat_z[i];
  }


  // likelihood
  bmd_z ~ normal(mu, sigma);
}
