#1: BMD~all
#gives direct effects of all predictors

#loop through calcification markers
#analysis with all calcification parameters, all inflammation parameters, all covariates and no sexes
calc_names <- c("all", colnames(dat$calc_markers))
fits_nosex <- list()

for (m in 1:length(calc_names)) {
  
  dat_nosex <- make_dat(
    data_to_use,
    imput = TRUE,
    sex = F,
    inflam_vars = "all",
    calc_vars = calc_names[m]
  )
  
  fits_nosex[[ calc_names[m] ]] <- cstan( file= "models/1bmd-all_i_c_2beta_imput.stan" , data=dat_nosex , chains=4, cores = 4, warmup = warmup, iter = iter, save_cmdstan_config=TRUE )
}


#adding sex differences
fits_c <- list()

for (m in 1:length(calc_names)) {
  
  dat_c <- make_dat(
    data_to_use,
    imput = TRUE,
    inflam_vars = "all",
    calc_vars = calc_names[m]
  )
  
  fits_c[[ calc_names[m] ]] <- cstan( file= "models/1bmd-all_i_c_2beta_imput.stan" , data=dat_c , chains=4, cores = 4, warmup = warmup, iter = iter )
}


#loop through inflammation markers. 
#no sex differences
inflam_names <- c("all", colnames(dat$inflam_markers))
fits_nosex_i <- list()

for (m in 1:length(inflam_names)) {
  
  dat_nosex <- make_dat(
    data_to_use,
    imput = FALSE,
    sex = F,
    inflam_vars = inflam_names[m],
    calc_vars = "all"
  )
  
  fits_nosex_i[[ inflam_names[m] ]] <- cstan( file= "models/1bmd-all_i_c_2beta.stan" , data=dat_nosex , chains=4, cores = 4, warmup = warmup, iter = iter )
}

# #with sex differences
# fits_i <- list()
# 
# for (m in 1:length(inflam_names)) {
#   
#   dat_i <- make_dat(
#     data_to_use,
#     imput = FALSE,
#     inflam_vars = inflam_names[m]
#   )
#   
#   fits_i[[ inflam_names[m] ]] <- cstan( file= "models/1bmd-all_i_c_2beta.stan" , data=dat_i , chains=4, cores = 4, warmup = warmup, iter = iter )
# }


#analysis with all calcification parameters, all inflammation parameters, all covariates, without imputing missing inflamation data from normal(0,1)
# dat_nosex_noimput <- make_dat(
#   data_to_use,
#   imput = FALSE,
#   sex = F
# )
# 
# fits_nosex_noimput <- cstan( file= "models/1bmd-all_i_c_2beta.stan" , data=dat_nosex_noimput , chains=4, cores = 4, warmup = 500, iter = 1000 )
# precis(fits_nosex_noimput, 3, c("alpha", "beta1", "beta2", "gamma", "epsilon", "eta", "theta", "sigma"))


dat_noimput <- make_dat(
  data_to_use,
  imput = FALSE
)

fits_noimput <- cstan( file= "models/1bmd-all_i_c_2beta.stan" , data=dat_noimput , chains=4, cores = 4, , warmup = warmup, iter = iter )
# precis(dat_noimput, 3, c("alpha", "beta1", "beta2", "gamma", "epsilon", "eta", "theta", "sigma"))


#2

dat_nosex <- make_dat(
  data_to_use,
  imput = TRUE,
  inflam_vars = "all",
  calc_vars = "all",
  sex = F
)

fits_2a <- cstan( file= "models/2a_bmd-nocalc_imput.stan" , data=dat_nosex , chains=4, cores = 4, warmup = warmup, iter = iter )
#precis(fits_2a, 3, c("alpha", "gamma", "epsilon", "eta", "theta", "sigma"))

fits_2b <- cstan( file= "models/2b_bmd-noinf.stan" , data=dat_nosex , chains=4, cores = 4, warmup = warmup, iter = iter )
#precis(fits_2b, 3, c("alpha", "beta1", "beta2", "gamma", "epsilon", "eta", "sigma"))


#2 by inflammation marker
inflam_names <- colnames(dat$inflam_markers)
fits_2a_i <- list()
for (m in 1:length(inflam_names)) {
  
  dat_nosex <- make_dat(
    data_to_use,
    imput = FALSE,
    sex = F,
    inflam_vars = inflam_names[m],
    calc_vars = "all"
  )
  
  fits_2a_i[[ inflam_names[m] ]] <- cstan( file= "models/2a_bmd-nocalc.stan" , data=dat_nosex , chains=4, cores = 4, warmup = warmup, iter = iter )
}


#3
dat_nosex <- make_dat(
  data_to_use,
  imput = TRUE,
  sex = F,
  inflam_vars = "all",
  calc_vars = "all"
)

fits_3 <- cstan( file= "models/3bmd_noinf_nocalc.stan" , data=dat_nosex , chains=4, cores = 4, warmup = warmup, iter = iter )
#precis(fits_3, 3, c("alpha", "gamma", "epsilon", "eta", "sigma"))


################
#calcification as outcome
###############

#4 & 5
# fits_4_nbz <- list()
# fits_4_pz <- list()
fits_4_pzln <- list()
fits_5 <- list()
calc_marks <- colnames(dat$calc_markers)
for (m in 1:length(calc_marks)) {

  dat_raw <- make_dat(
    data_to_use,
    imput = TRUE,
    inflam_vars = "all",
    calc_vars = "raw_calc",
    sex = F
  )

  dat_raw$calc <- dat_raw$calc_markers[,m]

  # fits_4_nbz[[ calc_marks[m] ]] <- cstan( file= "models/4calc-all_imput.stan" , data=dat_raw , chains=4, cores = 4, warmup = warmup, iter = iter )
  # fits_4_pz[[ calc_marks[m] ]] <- cstan( file= "models/4calc-all_imput_poisson.stan" , data=dat_raw , chains=4, cores = 4, warmup = warmup, iter = iter )
  fits_4_pzln[[ calc_marks[m] ]] <- cstan( file= "models/4_calc-all_zip_lognormal.stan" , data=dat_raw , chains=4, cores = 4, warmup = warmup, iter = iter )
  
  fits_5[[ calc_marks[m] ]] <- cstan( file= "models/5calc-noinf_poisson.stan" , data=dat_raw , chains=4, cores = 4, warmup = warmup, iter = iter )
}

# varying inflammation marker
fits_4_pzln_i <- list()
for (m in 1:length(calc_marks)) {
  for(i in 1:length(inflam_names)){
    dat_raw <- make_dat(
      data_to_use,
      imput = FALSE,
      inflam_vars = inflam_names[i],
      calc_vars = "raw_calc",
      sex = F
    )
    
    dat_raw$calc <- dat_raw$calc_markers[,m]
    
    fits_4_pzln_i[[ paste(calc_marks[m], "-", inflam_names[i]) ]] <- cstan( file= "models/4_calc-all_zip_lognormal.stan" , data=dat_raw , chains=4, cores = 4, warmup = warmup, iter = iter )
    

  }
}

# for (m in 1:length(calc_marks)) {
#   print(calc_marks[m])
#   print(precis(fits_4_nbz[[ calc_marks[m] ]], 2, c("alpha", "gamma", "epsilon", "eta", "theta", "alpha_zero", "gamma_zero", "epsilon_zero", "eta_zero", "theta_zero", "phi")))
#   print(precis(fits_4_pz[[ calc_marks[m] ]], 2, c("alpha", "gamma", "epsilon", "eta", "theta", "alpha_zero", "gamma_zero", "epsilon_zero", "eta_zero", "theta_zero")))
#   print(precis(fits_4_pzln[[ calc_marks[m] ]], 2, c("alpha", "gamma", "epsilon", "eta", "theta", "alpha_zero", "gamma_zero", "epsilon_zero", "eta_zero", "theta_zero")))
#   print(precis(fits_5[[ calc_marks[m] ]], 2, c("alpha", "gamma", "epsilon", "eta",  "alpha_zero", "gamma_zero", "epsilon_zero", "eta_zero")))
# }

# #6

inflam_names <- colnames(dat$inflam_markers)
fits_6 <- list()
for (m in 1:length(inflam_names)) {
  dat_i <- make_dat(
    data_to_use,
    imput = FALSE,
    inflam_vars = inflam_names[m],
    sex = F
  )
  dat_i$inflam_markers <- as.vector(dat_i$inflam_markers)

  fits_6[[ inflam_names[m] ]] <- cstan( file= "models/6infl-all.stan" , data=dat_i , chains=4, cores = 4, warmup = warmup, iter = iter )

}
# for (m in 1:length(inflam_names)) {
#   print(inflam_names[m])
#   print(precis(fits_6[[ inflam_names[m] ]], 2, c("alpha", "gamma", "epsilon", "eta", "sigma")))
# }
