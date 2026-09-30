data {
  
 // sample size
  int<lower=0> N;//n individuals
  int<lower=0> I;// n iflam markers
  int<lower=0> S;//one or two depending on whether calcification and age effects are sex specific

  // outcome
  vector[N] bmd_z;

  // predictors
  array[N] int<lower=1, upper=2> sex;
  vector[N] age_z;
  vector[N] ffm_z;
  vector[N] thoracic_fat_z;
  
  // inflammation indicators
  matrix[N, I] inflam_markers; // 1= crp, 2 = il1b, 3 = il6, 4 = tnfa

  }

parameters {
  // intercept
  vector[S] alpha;

  // regression coefficients
  vector[S] gamma;
  real epsilon;
  real eta;
  real theta;

  
  // residual SD
  real<lower=0> sigma;
  
}

transformed parameters {

  vector[N] inflam;
 
  inflam = inflam_markers * rep_vector(1.0 / I, I);

}

model {

  // priors
  alpha ~ normal(0,1);
  gamma ~ normal(0,1);
  epsilon ~ normal(0,1);
  eta ~ normal(0,1);
  theta ~ normal(0,1);
  sigma ~ exponential(1);
  
  
  // linear predictor
  vector[N] mu;

  for (i in 1:N) {

    mu[i] =
        alpha[sex[i]]
      + gamma[sex[i]] * age_z[i]
      + epsilon * ffm_z[i]
      + eta * thoracic_fat_z[i]
      + theta * inflam[i];
  }


  // likelihood
  bmd_z ~ normal(mu, sigma);
}
