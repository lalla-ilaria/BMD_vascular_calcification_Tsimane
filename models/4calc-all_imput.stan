data {

  int<lower=1> N;
  int<lower=0> I;// n iflam markers
  int<lower=1> S;

  array[N] int<lower=1,upper=S> sex;

  vector[N] age_z;
  vector[N] ffm_z;
  vector[N] thoracic_fat_z;
  // inflammation indicators
  matrix[N, I] inflam_markers; // 1= crp, 2 = il1b, 3 = il6, 4 = tnfa

  //missing inflammation
  int<lower=0> n_miss_inflam;
  array[n_miss_inflam] int miss_inflam;

  array[N] int<lower=0> calc;
}

parameters {

  // intercept
  vector[S] alpha;

  // regression coefficients
  vector[S] gamma;
  real epsilon;
  real eta;
  real theta;

  real<lower=0> phi;
  
  //zero inflated
  vector[S] alpha_zero;
  vector[S] gamma_zero;
  real epsilon_zero;
  real eta_zero;
  real theta_zero;

  // inflammation (Formative model)
  simplex[I] zeta;

  // Missing data
  vector[n_miss_inflam] il1b_miss;
  vector[n_miss_inflam]  il6_miss;
  vector[n_miss_inflam] tnfa_miss;
}

transformed parameters {

  vector[N] inflam;
  matrix[N,I] inflam_complete;

  inflam_complete = inflam_markers;

  for(j in 1:n_miss_inflam){
    int i = miss_inflam[j];

    inflam_complete[i,2] = il1b_miss[j];
    inflam_complete[i,3] = il6_miss[j];
    inflam_complete[i,4] = tnfa_miss[j];
  }
  
  inflam = inflam_complete * zeta;

}

model {

  // priors
  alpha ~ normal(0,1);
  gamma ~ normal(0,1);
  epsilon ~ normal(0,1);
  eta ~ normal(0,1);
  theta ~ normal(0,1);

  //zero inflated
  alpha_zero ~ normal(0,1);
  gamma_zero ~ normal(0,1);
  epsilon_zero ~ normal(0,1);
  eta_zero ~ normal(0,1);
  theta_zero ~ normal(0,1);
  
  //missing data
  il1b_miss ~ normal(0,1);
  il6_miss  ~ normal(0,1);
  tnfa_miss ~ normal(0,1);

  phi ~ exponential(1);

  for(i in 1:N){
    real log_mu;
    real logit_pi;
    real pi_i;

    log_mu =

        alpha[sex[i]]
      + gamma[sex[i]] * age_z[i]
      + epsilon * ffm_z[i]
      + eta * thoracic_fat_z[i]
      + theta * inflam[i];
      
    logit_pi =

        alpha_zero[sex[i]]
      + gamma_zero[sex[i]] * age_z[i]
      + epsilon_zero * ffm_z[i]
      + eta_zero * thoracic_fat_z[i]
      + theta_zero * inflam[i];
      

    pi_i = inv_logit(logit_pi);

    if(calc[i] == 0){

      target += log_sum_exp(
        bernoulli_lpmf(1 | pi_i),
        bernoulli_lpmf(0 | pi_i)
        + neg_binomial_2_log_lpmf(0 | log_mu, phi));
    } else {
      target +=
        bernoulli_lpmf(0 | pi_i)
        + neg_binomial_2_log_lpmf(
            calc[i] |
            log_mu,
            phi);
    }//if calc, else
  }//i in N

}

