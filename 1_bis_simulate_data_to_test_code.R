#simulate tsimane bmd

N <- 100
age_z <- rnorm(N, 0, 1)
sex <- rbinom(N, 1, 0.5)
diet <- rnorm(N, 0, 1)
ffm_z <- diet * 0.5 + age_z * - 0.3 + rnorm(N, 0, 0.3)
thoracic_fat_z <- diet * -0.5 + age_z * 0.3 + rnorm(N, 0, 0.3)
inflammation <- thoracic_fat_z * - 0.1 + rnorm(N, 0, 0.3)
log_crp_mean_z <- inflammation + rnorm(N, 0, 0.1)
log_il1b_cytokines_mean_z <- inflammation + rnorm(N, 0, 0.1)
log_il6_cytokines_mean_z <- inflammation + rnorm(N, 0, 0.1)
log_tnfa_cytokines_mean_z <- inflammation + rnorm(N, 0, 0.1)
calcification <- age_z * ifelse(sex, 0.3, 0.5) + inflammation * 0.1 + rnorm(N, 0, 0.3)
log_cac_z <- ifelse(calcification < 0, 0, calcification + rnorm(N, 0, 0.1))
log_tac_z <- calcification + rnorm(N, 0, 0.1)
log_aortic_arch_ca_z <- calcification + rnorm(N, 0, 0.1)
log_cac <- ifelse(calcification < 0, 0, calcification + rexp(N, 3))
log_tac <- ifelse(calcification < 0, 0, calcification + rexp(N, 3))
log_aac <- ifelse(calcification < 0, 0, calcification + rexp(N, 3))
cac <- (exp(log_cac)-1) * 20
tac <- (exp(log_tac)-1) * 60
aac <- (exp(log_aac)-1) * 100
bmd_z <- age_z * - ifelse(sex, 0.5, 0.3) + inflammation * - 0.1 + calcification * - 0.1 + rnorm(N, 0, 0.3)

dat_sim <- data.frame(
  age_z = age_z,
  sex = sex,
  ffm_z = ffm_z,
  thoracic_fat_z = thoracic_fat_z,
  bmd_z = bmd_z,
  log_crp_mean_z = log_crp_mean_z,
  log_il1b_cytokines_mean_z = log_il1b_cytokines_mean_z,
  log_il6_cytokines_mean_z = log_il6_cytokines_mean_z,
  log_tnfa_cytokines_mean_z = log_tnfa_cytokines_mean_z,
  log_cac = log_cac,
  log_tac = log_tac,
  log_aac = log_aac,
  cac = cac,
  tac = tac,
  aac = aac
)

dat <- make_dat(
  dat_sim, 
  imput = TRUE
)