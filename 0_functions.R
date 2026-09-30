##################################
# Function to make stan data list#
##################################
#prepare data for stan choosing which inflammation and calcification variables to consider, and divide sex
make_dat <- function(df,
                     imput = F,
                     sex = T,
                     inflam_vars = "all",
                     calc_vars = "all") {
  
  #removes missing data if not imputing
  if(!imput) {
    df <- df[complete.cases(df), ]
  } else {
    df <- df[complete.cases(df[,c("ffm_z", "thoracic_fat_z", "log_cac", "log_tac", "log_aac", "log_crp_mean_z")]), ]
  }
  
  miss_inflam <- which(
    is.na(df$log_il1b_cytokines_mean_z) |
      is.na(df$log_il6_cytokines_mean_z)  |
      is.na(df$log_tnfa_cytokines_mean_z)
  )
  
  #if imputing values, substitute NAs with 0, it will be overwritten in the model by normal(0,1), and remove data points missing other variables
  if(imput) {
    df$log_il1b_cytokines_mean_z[miss_inflam] <- 0
    df$log_il6_cytokines_mean_z[miss_inflam] <- 0
    df$log_tnfa_cytokines_mean_z[miss_inflam] <- 0
    
  }
  
  dat <- list (
    N = nrow(df),
    S = 2,
    sex = df$sex+1,
    bmd_z = df$bmd_z,
    age_z = df$age_z,
    ffm_z = df$ffm_z,
    thoracic_fat_z = df$thoracic_fat_z,
    miss_inflam = miss_inflam,
    n_miss_inflam = length(miss_inflam),
    inflam_markers = as.matrix(data.frame( crp = df$log_crp_mean_z,
                                           il1b = df$log_il1b_cytokines_mean_z,
                                           il6 = df$log_il6_cytokines_mean_z,
                                           tnfa = df$log_tnfa_cytokines_mean_z)),
    #crp = df$log_crp_mean_z,
    #il1b = df$log_il1b_cytokines_mean_z,
    #il6 = df$log_il6_cytokines_mean_z,
    #tnfa = df$log_tnfa_cytokines_mean_z,
    calc_markers = as.matrix(data.frame( cac = df$log_st_cac,
                                         tac = df$log_st_tac,
                                         aac = df$log_st_aac)),
    any_calc = as.integer(df$log_cac > 0 |
                            df$log_tac > 0 |
                            df$log_aac > 0)
  )
  #substitutes standardized calcifications with raw values for certain analyses
  if(calc_vars == "raw_calc") dat$calc_markers = as.matrix(data.frame( cac = df$cac,
                                                                       tac = df$tac,
                                                                       aac = df$aac))
  
  
  if(sex == F) {
    dat$S <- 1
    dat$sex <- rep(1, dat$N)
  }
  
  # inflammation markers
  if (inflam_vars %in% c("crp", "il1b", "il6", "tnfa")) {
    
    dat$inflam_markers <-
      dat$inflam_markers[, inflam_vars, drop = FALSE]
  }
  
  # calcification markers
  if (calc_vars %in% c("cac", "tac", "aac")) {
    
    dat$calc_markers <-
      dat$calc_markers[, calc_vars, drop = FALSE]
    dat$any_calc <-
      dat$calc_markers[, calc_vars] > 0 #any calcification on the marker of interest, not on any marker
  }
  
  # dimensions for Stan
  dat$I <- ncol(dat$inflam_markers)
  dat$C <- ncol(dat$calc_markers)
  
  #make sure correct types of data for Stan
  dat$N <- as.integer(dat$N)
  dat$I <- as.integer(dat$I)
  dat$C <- as.integer(dat$C)
  dat$S <- as.integer(dat$S)
  
  dat$sex <- as.integer(dat$sex)
  dat$miss_inflam <- as.integer(dat$miss_inflam)
  dat$n_miss_inflam <- as.integer(dat$n_miss_inflam)
  
  dat$any_calc <- as.integer(dat$any_calc)
  
  dat$bmd_z <- as.numeric(dat$bmd_z)
  dat$age_z <- as.numeric(dat$age_z)
  dat$ffm_z <- as.numeric(dat$ffm_z)
  dat$thoracic_fat_z <- as.numeric(dat$thoracic_fat_z)
  dat$inflam_markers <- as.matrix(dat$inflam_markers)
  dat$calc_markers <- as.matrix(dat$calc_markers)
  return(dat)
}

################################################
# Functions to prep SI table on diff btw groups#
################################################
mean_sd <- function(x){
  
  sprintf(
    "%.1f (%.1f)",
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE)
  )
}

percent_yes <- function(x){
  
  sprintf(
    "%.0f",
    100 * mean(x, na.rm = TRUE)
  )
}

make_row <- function(dat, var, label, binary = FALSE){
  
  g1 <- subset(dat, group == "low bmd any calc mild")
  g2 <- subset(dat, group == "low bmd any calc mod.sev")
  g3 <- subset(dat, group == "other")
  
  if(binary){
    
    s1 <- percent_yes(g1[[var]])
    s2 <- percent_yes(g2[[var]])
    s3 <- percent_yes(g3[[var]])
    
    m1 <- 100 * mean(g1[[var]], na.rm = TRUE)
    m2 <- 100 * mean(g2[[var]], na.rm = TRUE)
    m3 <- 100 * mean(g3[[var]], na.rm = TRUE)
    
  } else {
    
    s1 <- mean_sd(g1[[var]])
    s2 <- mean_sd(g2[[var]])
    s3 <- mean_sd(g3[[var]])
    
    m1 <- mean(g1[[var]], na.rm = TRUE)
    m2 <- mean(g2[[var]], na.rm = TRUE)
    m3 <- mean(g3[[var]], na.rm = TRUE)
    
  }
  
  data.frame(
    Variable = label,
    Group1 = s1,
    Group2 = s2,
    Group3 = s3,
    `1 vs 2` = round(abs(m1 - m2),1),
    `1 vs 3` = round(abs(m1 - m3),1),
    `2 vs 3` = round(abs(m2 - m3),1)
  )
}
####################################
# Function to print table of precis#
####################################
#rename parameters
rename_pars <- function(x){
  
  x <- gsub("alpha\\[1\\]", "α (Intercept,sex 1)", x)
  x <- gsub("alpha\\[2\\]", "α (Intercept,sex 2)", x)
  
  x <- gsub("beta1\\[1\\]", "β₁ (Calcification presence, sex 1)", x)
  x <- gsub("beta1\\[2\\]", "β₁ (Calcification presence, sex 2)", x)
  
  x <- gsub("beta2\\[1\\]", "β₂ (Calcification amount, sex 1)", x)
  x <- gsub("beta2\\[2\\]", "β₂ (Calcification amount, sex 2)", x)
  
  x <- gsub("gamma\\[1\\]", "γ (Age, sex 1)", x)
  x <- gsub("gamma\\[2\\]", "γ (Age, sex 2)", x)
  x <- gsub("epsilon$", "ε (FFM)", x)
  x <- gsub("theta$", "θ (Inflammation)", x)
  x <- gsub("eta$", "η (Thoracic fat)", x)
  
  # x <- gsub("xi\\[1\\]", "ξ (CAC weight)", x)
  # x <- gsub("xi\\[2\\]", "ξ (TAC weight)", x)
  # x <- gsub("xi\\[3\\]", "ξ (AAC weight)", x)
  # 
  # x <- gsub("zeta\\[1\\]", "ζ (CRP weight)", x)
  # x <- gsub("zeta\\[2\\]", "ζ (IL-1β weight)", x)
  # x <- gsub("zeta\\[3\\]", "ζ (IL6 weight)", x)
  # x <- gsub("zeta\\[4\\]", "ζ (TNF-α weight)", x)
  x <- gsub("alpha_zero\\[1\\]", "α₀ (ZI intercept,women)", x)
  x <- gsub("alpha_zero\\[2\\]", "α₀ (ZI intercept, men)", x)
  x <- gsub("gamma_zero\\[1\\]", "γ₀ (ZI Age, women)", x)
  x <- gsub("gamma_zero\\[2\\]", "γ₀ (ZI Age, men)", x)
  x <- gsub("epsilon_zero$", "ε₀ (ZI FFM)", x)
  x <- gsub("theta_zero$", "θ₀ (ZI Inflammation)", x)
  x <- gsub("eta_zero$", "η₀ (ZI Thoracic fat)", x)
  x <- gsub("phi$", "φ (Overdispersion)", x)
  
  x
}

extract_precis <- function(fit, pars){
  
  p <- precis(fit, depth = 2, pars = pars)
  
  out <- data.frame(
    parameter = rownames(p),
    mean = p[, "mean"],
    lower = p[, "5.5%"],
    upper = p[, "94.5%"],
    row.names = NULL
  )
  
  sig <- (out$lower > 0 & out$upper > 0) |
    (out$lower < 0 & out$upper < 0)
  
  est <- sprintf(
    "%.2f [%.2f, %.2f]",
    out$mean, out$lower, out$upper
  )
  
  #est[sig] <- paste0("**", est[sig], "**")
  
  data.frame(
    parameter = rename_pars(out$parameter),
    estimate = est
  )
}


##################################
# Functions to compare parameters#
##################################
pretty_par_name <- function(par){
  
  switch(
    par,
    alpha = "α (Intercept)",
    beta1 = "β (Calcification present)",
    beta2 = "β (Calcification amount)",
    gamma = "γ (Age)",
    epsilon = "ε (FFM)",
    eta = "η (Thoracic fat)",
    theta = "θ (Inflammation composite)",
    phi = "φ (Overdispersion)",
    # zeta = "ζ (Inflammation weights)",
    # xi = "ξ (Calcification weights)",
    alpha_zero = "α₀ (ZI intercept)",
    gamma_zero = "γ₀ (ZI age)",
    epsilon_zero = "ε₀ (ZI FFM)",
    eta_zero = "η₀ (ZI thoracic fat)",
    theta_zero = "θ₀ (ZI inflammation)",
    par
  )
}

#shows direct comparison btw parameters from two fitted posteriors
compare_param <- function(post1, post2, par, post1_label = "First model", post2_label = "Second model") {
  
  d1 <- as.vector(post1[[par]])
  d2 <- as.vector(post2[[par]])
  
  data.frame(
    model = c(post1_label, post2_label),
    mean  = c(mean(d1), mean(d2)),
    sd    = c(sd(d1), sd(d2)),
    lower = c(PI(d1)[1], PI(d2)[1]),
    upper = c(PI(d1)[2], PI(d2)[2])
  )
}

#plots densities parameters
compare_density <- function(
    post1,
    post2,
    par,
    cols = c(femalecol, malecol),
    ltys = c(1, 2)
){
  label <- pretty_par_name(par)#renames to nice name
  
  x1 <- post1[[par]]
  x2 <- post2[[par]]
  if( dim(x1)[2] != dim(x2)[2]) warning("Models are fitted with different number of sexes")
  # matrix parameter (sex-specific)
  if(dim(x1)[2] == 2){
    ncol_par <- ncol(x1)
    
    for(j in 1:ncol_par){
      
      d1 <- density(x1[,j])
      d2 <- density(x2[,j])
      
      if(j == 1){
        plot(d1, lwd = 2, col = cols[j],
             main = label, xlab = label, xlim = c(-1,1), ylim = c(0,9))
      } else {
        lines(d1, lwd = 2, col = cols[j])
      }
      lines(d2, lwd = 2, col = cols[j], lty = ltys[2])
    }
    
    legend("topright",
           legend = c("Women complete", "Men complete",
             "Women imputed", "Men imputed"),
           col = rep(cols[1:ncol_par], 2),
           lty = c(rep(1, ncol_par), rep(2, ncol_par)),
           lwd = 2, bty = "n" )
    
  } else { #if single sex
    
    d1 <- density(as.vector(x1))
    d2 <- density(as.vector(x2))
    
    plot(d1, lwd = 2, main = label, xlab = label, xlim = c(-1,1), ylim = c(0,9), col = allcol)
    
    lines(d2, lwd = 2, lty = 2, col = allcol)
    
    legend("topright",
           legend = c("complete", "imputed"),
           lwd = 2, lty = c(1,2), bty = "n", col = allcol )
  }
}


############################################
#Function to calculate correlation matrixes#
############################################
panel.cor <- function(
    x,
    y,
    digits = 2,
    cex.cor = 1.2,
    method = "pearson",
    return_value = FALSE,
    ...
){
  
  r <- cor(
    x,
    y,
    use = "pairwise.complete.obs",
    method = method
  )
  
  if(return_value){
    return(round(r, digits))
  }
  
  usr <- par("usr")
  on.exit(par(usr))
  
  par(usr = c(0,1,0,1))
  
  text(
    0.5,
    0.5,
    paste0("r = ", round(r, digits)),
    cex = cex.cor
  )
}

##################################################
# Functions to build generated quantities entries#
##################################################

gq_entry <- function(interpretation, draws = NULL, mean = NULL, pi_width = 0.89, digits = 2, section = "ungrouped") {
  if(!is.null(draws)){
    m <- round(mean(draws), digits)
    pi <- round(PI(draws, prob = pi_width), digits)
  } else {
    m <- round(mean, digits)
    pi <- NULL
  }
  list(interpretation = interpretation, mean = m, PI = pi, section = section)
}

gq_add <- function(gq, name, ...) {
  if(name %in% names(gq)) warning("Overwriting existing gq entry: ", name)
  gq[[name]] <- gq_entry(...)
  gq
}
gq_export_json <- function(gq, path = "gq.json") {
  records <- lapply(names(gq), function(nm) {
    e <- gq[[nm]]
    list(
      name = nm,
      interpretation = e$interpretation,
      mean = e$mean,
      PI_lower = if(!is.null(e$PI)) e$PI[1] else NA,
      PI_upper = if(!is.null(e$PI)) e$PI[2] else NA,
      section = if(!is.null(e$section)) e$section else "ungrouped"
    )
  })
  jsonlite::write_json(records, path, auto_unbox = TRUE, na = "null")
  invisible(path)
}

######################################################
# Functions to produce a single counterfactual value #
######################################################
#calculates posterior distribution of BMD for counterfactual (average) individual (male) 
#note that here is the function that has to match the actual function used in the model
#current specification allows to change such function depending on the parameters present in the model fit given to the function

#BMD as outcome
counterfactual_BMD <- function(post, sex = 1,
                               calc = 0,
                               age = counterfactual_age_zscore,
                               ffm = 0,
                               thoracic_fat = 0,
                               inflam = c(0, 0, 0, 0),
                               calc_composite = FALSE,
                               inflam_composite = FALSE
){
  
  n_draws <- nrow(post$alpha)
  #calculates composites
  if("theta" %in% names(post)){
    if(inflam_composite){
      c_i <- inflam #run over composite inflammation
    } else {
      c_i <- rowMeans(
        matrix(
          rep(inflam, each = n_draws),
          nrow = n_draws,
          ncol = length(inflam)))
    }#individual inflammation markers
  }#if model includes inflammation
  if("beta1" %in% names(post)){
    if(calc_composite){
      c_c <- calc
      presence <- 1
    } else { 
      c_c <- rowMeans(
        matrix(
          rep(calc, each = n_draws),
          nrow = n_draws,
          ncol = length(calc)))
      presence <- as.numeric(sum(calc) > 0)
    }#individual calcification markers
  }#if model includes calcification
  #full model for comparison
  # distribution <-
  #   post$alpha[, sex] +
  #   post$beta1[, sex] * presence +
  #   post$beta2[, sex] * c_c +
  #   post$gamma[, sex] * age +
  #   post$epsilon * ffm +
  #   post$eta * thoracic_fat +
  #   post$theta * c_i
  
  #flexible model that includes parameters present in posteriors passed to function
  distribution <- post$alpha[, sex] 
  
  if("beta1" %in% names(post))    distribution <- distribution + post$beta1[,sex] * presence + post$beta2[,sex] * c_c
  if("theta" %in% names(post))    distribution <- distribution + post$theta * c_i
  if("epsilon" %in% names(post))  distribution <- distribution + post$epsilon * ffm
  if("eta" %in% names(post))      distribution <- distribution + post$eta * thoracic_fat
  if("gamma" %in% names(post))    distribution <- distribution + post$gamma[,sex] * age
  
  return(distribution)
}

#calc as outcome
#note that model is negative binomial, so it returns expected calcification for 


counterfactual_calc <- function(
    post,
    sex = 1,
    age = counterfactual_age_zscore,
    ffm = 0,
    thoracic_fat = 0,
    inflam = c(0,0,0,0),
    inflam_composite = FALSE
){
  
  n_draws <- nrow(post$alpha)
  
  # inflammation composite
  if("theta" %in% names(post)){
    if(inflam_composite){
       c_i <- inflam
    } else {
      c_i <- rowMeans(
        matrix(
          rep(inflam, each = n_draws),
          nrow = n_draws,
          ncol = length(inflam))
      )}}
  
  # count process

  log_mu <- post$alpha[,sex]
  
  if("theta" %in% names(post))    log_mu <- log_mu + post$theta * c_i
  if("epsilon" %in% names(post))  log_mu <- log_mu + post$epsilon * ffm
  if("eta" %in% names(post))      log_mu <- log_mu + post$eta * thoracic_fat
  if("gamma" %in% names(post))    log_mu <- log_mu + post$gamma[,sex] * age
  
  if("sigma_u" %in% names(post)) {
    mu <- exp(log_mu + 0.5 * post$sigma_u^2)
  } else {
    mu <- exp(log_mu)
  }
  # zero process

  if("alpha_zero" %in% names(post)){
    
    logit_pi <- post$alpha_zero[,sex]
    
    if("theta_zero" %in% names(post))   logit_pi <- logit_pi + post$theta_zero * c_i
    if("epsilon_zero" %in% names(post)) logit_pi <- logit_pi + post$epsilon_zero * ffm
    if("eta_zero" %in% names(post))     logit_pi <- logit_pi + post$eta_zero * thoracic_fat
    if("gamma_zero" %in% names(post))   logit_pi <- logit_pi + post$gamma_zero[,sex] * age
    
    pi <- plogis(logit_pi)
    distribution <- (1 - pi) * mu
    
  } else {
    distribution <- mu
  }
  
  distribution
}


#infl as outcome
counterfactual_infl <- function(post, sex = 1,
                                age = counterfactual_age_zscore,
                                ffm = 0,
                                thoracic_fat = 0
){
  
  distribution <-
    post$alpha[, sex] +
    post$gamma[, sex] * age +
    post$epsilon * ffm +
    post$eta * thoracic_fat 
  
  return(distribution)
}

###############################################
# Function for calculating counterfactual line#
###############################################

counterfactual_line <- function(
    post,
    outcome = "BMD",
    vary,
    vary_seq,
    sex = 1, #sex = 0 means average between the sexes weighing by counts, 1 is female 2 is male
    covariates
) {
  
  counterfactual_fn <- switch(
    outcome,
    BMD  = counterfactual_BMD,
    calc = counterfactual_calc,
    infl = counterfactual_infl,
    stop("Unknown outcome")
  )
  out <- lapply(vary_seq, function(x) {
    
    # default values
    args <- covariates
    args$post <- post

    # replace selected variable
    args[[vary]] <- x
    
    # posterior predictions
    #averaging between sexes
    if(sex == 0){
      #assigns sex
      args_f <- args
      args_f$sex <- 1
      
      args_m <- args
      args_m$sex <- 2
      #calculates weights by sex
      w_f <- mean(dat$sex == 1)
      w_m <- mean(dat$sex == 2)
      #calculates counterfactual by sex
      mu_f <- do.call(counterfactual_fn, args_f)
      mu_m <- do.call(counterfactual_fn, args_m)
      #averages between the sexes weighing by counts
      mu <- w_f * mu_f + w_m * mu_m
      
    } else {
    #for sex as called by sex = 1|2  
      args$sex <- sex
      mu <- do.call(counterfactual_fn, args)
      
    }
    #creates data frame with mean and 89%PI
    data.frame(
      x = if(length(x) == 1) x else mean(x),
      mean = mean(mu),
      lower = PI(mu)[1],
      upper = PI(mu)[2]
    )
  })
  
  do.call(rbind, out)
}
##############################
# Show estimates for presence#
##############################

# Calculate jump from zero to one calcification
calc_jump <- function(post, sex = 1,
                      age = counterfactual_age_zscore,
                      ffm = 0,
                      thoracic_fat = 0,
                      inflam = 0){
  
  mu0 <- counterfactual_BMD(
    post = post,
    sex = sex,
    calc = 0,
    age = age,
    ffm = ffm,
    thoracic_fat = thoracic_fat,
    inflam = inflam
  )
  
  mu1 <- counterfactual_BMD(
    post = post,
    sex = sex,
    calc = 1,   # “any calcification”
    age = age,
    ffm = ffm,
    thoracic_fat = thoracic_fat,
    inflam = inflam
  )
  
  list(
    no_calc = mu0,
    any_calc = mu1,
    diff = mu1 - mu0
  )
}

################################
# Plots to show counterfactuals#
################################

plot_counterfactual_sex <- function(
    post,
    outcome = "BMD",
    vary,
    vary_seq = seqalong,
    covariates = NULL,
    xlim = NULL,
    ylim = NULL,
    femalecol = NULL,
    malecol = NULL,
    main = "",
    xlab = NULL,
    ylab = NULL,
    latent_points = FALSE,
    latent_var = NULL,
    dat = dat
){
  
  if (is.null(covariates)) covariates <- get("covariates", envir = .GlobalEnv)
  if (is.null(xlim))       xlim       <- get("xlim", envir = .GlobalEnv)
  if (is.null(ylim))       ylim       <- get("ylim", envir = .GlobalEnv)
  if (is.null(femalecol))  femalecol  <- get("femalecol", envir = .GlobalEnv)
  if (is.null(malecol))    malecol    <- get("malecol", envir = .GlobalEnv)
  if (is.null(ylab))       ylab       <- paste("Expected", outcome, "(z-scored)")
  if (is.null(xlab))       xlab       <- vary
  
  if(outcome == "calc") ylab <- "Expected calcification (log)"
  #add jump to side of calcification plots
  if(vary == "calc" & !latent_points) jump_space <- diff(xlim) * 0.23 else jump_space <- 0

  d_f <- counterfactual_line(
    post = post,
    outcome = outcome,
    sex = 1,
    covariates = covariates,
    vary = vary,
    vary_seq = vary_seq
  )
  
  d_m <- counterfactual_line(
    post = post,
    outcome = outcome,
    sex = 2,
    covariates = covariates,
    vary = vary,
    vary_seq = vary_seq
  )
  
  plot(
    d_f$x, d_f$mean,
    type = "n",
    xlim = c(xlim[1] - jump_space, xlim[2]), 
    ylim = ylim,
    xlab = xlab, ylab = ylab,
    main = main
  )
  
  abline(h = 0, col = "grey80", lwd = 1.5, lty = 3)
  abline(h = 100, col = "grey80", lwd = 1.5, lty = 3)
  abline(h = 400, col = "grey80", lwd = 1.5, lty = 3)
  text(max(xlim) - 1, 50, 'min risk')
  text(max(xlim) - 1, 200, 'low risk')
  text(max(xlim) - 1, 500, 'high risk')
  
  
  if(latent_points){
    latent_mean <- colMeans(post[[latent_var]])
    latent_low  <- apply(post[[latent_var]], 2, PI)[1,]
    latent_high <- apply(post[[latent_var]], 2, PI)[2,]
    for(i in seq_along(latent_mean)){
      segments(
        latent_low[i], dat$bmd_z[i], latent_high[i], dat$bmd_z[i],
        col = adjustcolor(cols[dat$sex[i]], 0.3)
      )
    }#latentpoints
    
    points(
      latent_mean, dat$bmd_z,
      pch = 16, col = adjustcolor(cols[dat$sex], 0.5))
  }
  
  polygon(
    c(d_f$x, rev(d_f$x)), c(d_f$lower, rev(d_f$upper)),
    border = NA, col = adjustcolor(femalecol, alpha.f = 0.2))
  
  polygon(
    c(d_m$x, rev(d_m$x)), c(d_m$lower, rev(d_m$upper)),
    border = NA, col = adjustcolor(malecol, alpha.f = 0.2)
  )
  
  lines(d_f$x, d_f$mean, lwd = 2, col = femalecol)
  lines(d_m$x, d_m$mean, lwd = 2, col = malecol)

    #add jump from zero to one calcification
  if(vary == "calc" & !latent_points){
    abline(v = 0, col = "grey80", lwd = 1.5, lty = 3)
    #abline(v = log10(101), col = "grey80", lwd = 1.5, lty = 4)
    #abline(v = log10(401), col = "grey80", lwd = 1.5, lty = 4)
    mu0_f  <- counterfactual_BMD(post, sex = 1, calc = covariates$calc, age = covariates$age, ffm = covariates$ffm, thoracic_fat = covariates$thoracic_fat, inflam = covariates$inflam) 
    mu0_m  <- counterfactual_BMD(post, sex = 2, calc = covariates$calc, age = covariates$age, ffm = covariates$ffm, thoracic_fat = covariates$thoracic_fat, inflam = covariates$inflam) 
    xjump_f <- xlim[1] - jump_space * 0.8
    xjump_m <- xlim[1] - jump_space * 0.4
    segments(xjump_f, PI(mu0_f)[1], xjump_f,PI(mu0_f)[2],
      lwd = 5, col = adjustcolor(femalecol, alpha.f = 0.2))
    
    points( xjump_f, mean(mu0_f),
      pch = 15, col = femalecol)
    
    segments(xjump_m, PI(mu0_m)[1],
      xjump_m, PI(mu0_m)[2],
      lwd = 5, col = adjustcolor(malecol, alpha.f = 0.2))
    
    points( xjump_m, mean(mu0_m),
      pch = 15, col = malecol)
    
    text(-0.42, 0.9, "No \n calcification")
    text(1, 0.9, "Calcification \n severity")
  }

  invisible(list(female = d_f, male = d_m))
}


#plot single sex function
plot_counterfactual <- function(
    post,
    outcome = "BMD",
    vary,
    vary_seq = seqalong,
    covariates = NULL,
    xlim = NULL,
    ylim = NULL,
    col = allcol,
    main = "",
    xlab = NULL,
    ylab = NULL,
    latent_points = FALSE,
    latent_var = NULL,
    dat = dat,
    add = FALSE
){
  
  if (is.null(covariates)) covariates <- get("covariates", envir = .GlobalEnv)
  if (is.null(xlim))       xlim       <- get("xlim", envir = .GlobalEnv)
  if (is.null(ylim))       ylim       <- get("ylim", envir = .GlobalEnv)
  if (is.null(col))        col        <- "black"
  if (is.null(ylab))       ylab       <- paste("Expected", outcome, "(z-scored)")
  if (is.null(xlab))       xlab       <- vary
  if(outcome == "calc") ylab <- "Expected calcification (log)"
  if(vary == "calc") jump_space <- diff(xlim) * 0.23 else jump_space <- 0
  
  d <- counterfactual_line(
    post = post,
    outcome = outcome,
    covariates = covariates,
    vary = vary,
    vary_seq = vary_seq
  )
  
  if(!add) {
    plot(
      d$x, d$mean, type = "n",
      xlim = c(xlim[1] - jump_space, xlim[2]), ylim = ylim,
      xlab = xlab, ylab = ylab, 
      main = main)
  }
  abline(h = 0, col = "grey80", lwd = 1.5, lty = 3)
  
  if(latent_points){
    latent_mean <- colMeans(post[[latent_var]])
    latent_low  <- apply(post[[latent_var]], 2, PI)[1,]
    latent_high <- apply(post[[latent_var]], 2, PI)[2,]
    for(i in seq_along(latent_mean)){
      segments(
        latent_low[i], dat$bmd_z[i],latent_high[i], dat$bmd_z[i],
        col = adjustcolor(col, 0.3))
    }
    points(
      latent_mean, dat$bmd_z,
      pch = 16, col = adjustcolor(col, 0.5))
  }#latent points
  
  polygon(
    c(d$x, rev(d$x)), c(d$lower, rev(d$upper)),
    border = NA, col = adjustcolor(col, alpha.f = 0.2))
  
  lines(
    d$x,d$mean,
    lwd = 2, col = col)

  if(vary == "calc"){
    abline(v = 0, col = "grey80", lwd = 1.5, lty = 3)
    #abline(v = log10(101), col = "grey80", lwd = 1.5, lty = 4)
    #abline(v = log10(401), col = "grey80", lwd = 1.5, lty = 4)
    mu0  <- counterfactual_BMD(post, sex = 1, calc = covariates$calc, age = covariates$age, ffm = covariates$ffm, thoracic_fat = covariates$thoracic_fat, inflam = covariates$inflam) 
    xjump <- xlim[1] - jump_space + 0.35
    segments(xjump, PI(mu0)[1], xjump,PI(mu0)[2],
             lwd = 5, col = adjustcolor(allcol, alpha.f = 0.2))
    
    points( xjump, mean(mu0),
            pch = 15, col = allcol)
    
    text(-0.42, 0.9, "No \n calcification")
    text(1, 0.9, "Calcification \n severity")
  }
  
  invisible(d)
}


#legend helper
legend_ribbon <- function(
    x0, y0,
    label,
    col,
    alpha = 0.2,
    lty = 1,
    width = 0.4,
    height = 0.08
){
  rect(
    x0, y0, x0 + width, y0 + height,
    col = adjustcolor(col, alpha.f = alpha), border = NA)
  segments(
    x0, y0 + height/2, x0 + width, y0 + height/2,
    col = col, lwd = 2, lty = lty )
  text(
    x0 + width + 0.05, y0 + height/2,
    label,adj = 0)
}

###########################################################
# Posterior prediction checks for fit to negative binomial#
###########################################################
ppd_calc <- function(post, dat){
  n_draws <- length(post$phi)
  N <- dat$N
  yrep <- matrix(NA, nrow = n_draws, ncol = N)
  
  for(d in seq_len(n_draws)){
    # count process
    log_mu <- post$alpha[d] + post$gamma[d] * dat$age_z + post$epsilon[d] * dat$ffm_z +
      post$eta[d] * dat$thoracic_fat_z + post$theta[d] * post$inflam[d,]
    mu <- exp(log_mu)
    
    # zero process
    if("alpha_pi" %in% names(post)){
      logit_pi <- post$alpha_pi[d] + post$gamma_pi[d] * dat$age_z + post$epsilon_pi[d] * dat$ffm_z +
        post$eta_pi[d] * dat$thoracic_fat_z + post$theta_pi[d] * post$inflam[d,]
      pi <- plogis(logit_pi)
    } else {
      pi <- rep(0, N)
    }
    
    # simulate
    structural_zero <- rbinom(N, size = 1, prob = pi)
    nb_draw <- rnbinom(N, mu = mu, size = post$phi[d])
    yrep[d,] <- ifelse(structural_zero == 1, 0, nb_draw)
  }
  yrep
}

#for distributions
plot_ppc_distribution <- function(y, yrep){
  hist(y, breaks = 50, freq = FALSE, main = "", xlab = "Calcification")
  for(i in sample(nrow(yrep), 50)) lines(density(yrep[i,]), col = adjustcolor(partialcol1, 0.3))
  lines(density(y), lwd = 3, col = totalcol)
}

#for proportions of zero
ppc_zeros <- function(y, yrep){
  obs_zero <- mean(y == 0)
  rep_zero <- apply(yrep, 1, function(x) mean(x == 0))
  hist(rep_zero, breaks = 30, main = "Proportion of zeros", xlab = "Simulated", 
       col = partialcol1, border = partialcol1)
  abline(v = obs_zero, lwd = 3, col = totalcol)
}

#to compare mean and variance
ppc_meanvar <- function(y, yrep){
  rep_mean <- rowMeans(yrep)
  rep_var <- apply(yrep, 1, var)
  plot(log10(rep_mean + 1), log10(rep_var + 1), pch = 16, 
       col = adjustcolor(partialcol2, 0.4),
       xlab = "Mean", ylab = "Variance",
       main = paste("Real data mean:", round(log10(mean(y) + 1), 2), "var:", round(log10(var(y) + 1), 2)))
  points(log10(mean(y) + 1), log10(var(y) + 1), pch = 16, cex = 2)
}

###########################################################
# Posterior prediction checks for fit to zero-inflated Poisson
###########################################################
ppd_calc_zip <- function(post, dat){
  n_draws <- nrow(post$alpha)
  N <- dat$N
  yrep <- matrix(NA, nrow = n_draws, ncol = N)
  
  for(d in seq_len(n_draws)){
    log_mu <- post$alpha[d] + post$gamma[d] * dat$age_z + post$epsilon[d] * dat$ffm_z +
      post$eta[d] * dat$thoracic_fat_z + post$theta[d] * post$inflam[d,]
    mu <- exp(log_mu)
    
    logit_pi <- post$alpha_zero[d] + post$gamma_zero[d] * dat$age_z + post$epsilon_zero[d] * dat$ffm_z +
      post$eta_zero[d] * dat$thoracic_fat_z + post$theta_zero[d] * post$inflam[d,]
    pi <- plogis(logit_pi)
    
    structural_zero <- rbinom(N, size = 1, prob = pi)
    poisson_draw <- rpois(N, lambda = mu)
    yrep[d,] <- ifelse(structural_zero == 1, 0, poisson_draw)
  }
  yrep
}

###########################################################
# Posterior prediction checks for fit to zero-inflated lognormal-Poisson mixture
###########################################################
ppd_calc_zipln <- function(post, dat){
  n_draws <- nrow(post$alpha)
  N <- dat$N
  yrep <- matrix(NA, nrow = n_draws, ncol = N)
  
  for(d in seq_len(n_draws)){
    log_mu <- post$alpha[d] + post$gamma[d] * dat$age_z + post$epsilon[d] * dat$ffm_z +
      post$eta[d] * dat$thoracic_fat_z + post$theta[d] * post$inflam[d,]
    
    # draw the latent overdispersion term and add it, then simulate Poisson
    u <- rnorm(N, mean = 0, sd = post$sigma_u[d])
    mu <- exp(log_mu + u)
    
    logit_pi <- post$alpha_zero[d] + post$gamma_zero[d] * dat$age_z + post$epsilon_zero[d] * dat$ffm_z +
      post$eta_zero[d] * dat$thoracic_fat_z + post$theta_zero[d] * post$inflam[d,]
    pi <- plogis(logit_pi)
    
    structural_zero <- rbinom(N, size = 1, prob = pi)
    poisson_draw <- rpois(N, lambda = mu)
    yrep[d,] <- ifelse(structural_zero == 1, 0, poisson_draw)
  }
  yrep
}
################################################################
# Function to prepare values for tables by inflammation markers# 
################################################################
build_infl_tab <- function(marker, pars = c("alpha","gamma","epsilon","eta","theta")) {
  full   <- extract_precis(fits_nosex_i[[marker]], pars) |> rename("Full model" = estimate)
  nocalc <- extract_precis(fits_2a_i[[marker]],    pars) |> rename("No calc model" = estimate)
  full_join(nocalc, full, by = "parameter")
}
make_inflam_table <- function(inflam){
  
  tab <- NULL
  
  for(v in calc_vars){
    
    tmp <-
      extract_precis(
        fits_4_pzln_i[[paste(v, "-", inflam)]],
        pars
      ) |>
      rename(
        !!paste("Effect on", toupper(v)) := estimate
      )
    
    if(is.null(tab)){
      tab <- tmp
    } else {
      tab <- full_join(tab, tmp, by = "parameter")
    }
  }
  
  tab
}