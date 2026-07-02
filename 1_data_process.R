data_expand <- read.csv("bmd_data_rev19_r.csv")
# Define as factors groups used in Table S28
data_expand$group_lowbmd_cac_level <- as.factor(data_expand$group_lowbmd_cac_level)
data_expand$group_lowbmd_tac_level <- as.factor(data_expand$group_lowbmd_tac_level)
data_expand$group_lowbmd_aac_level <- as.factor(data_expand$group_lowbmd_aac_level)
data_expand$group_lowbmd_any_calc_level <- as.factor(data_expand$group_lowbmd_any_calc_level)

# For exploring potential quadratic effects
data_expand$age2 = data_expand$age * data_expand$age
data_expand$n_births2 = data_expand$n_births * data_expand$n_births

# Create diverse categories
# data_expand$bmi_cat <- cut(
#   data_expand$bmi,
#   breaks = c(-Inf, 18.5, 25, 30, Inf),
#   labels = c("Underweight", "Normal weight", "Overweight", "Obesity"),
#   right = FALSE
# )

#data_expand$bmi_cat = as.factor(data_expand$bmi_cat)

# Z-score variables: pooled sexes
data_expand$bmd_z <- scale(data_expand$bmd, center = TRUE, scale = TRUE)
data_expand$age_z <- scale(data_expand$age, center = TRUE, scale = TRUE)
data_expand$ffm_z <- scale(data_expand$fatfreemass, center = TRUE, scale = TRUE)
data_expand$thoracic_fat_z <- scale(data_expand$thoracic_fat, center = TRUE, scale = TRUE)
data_expand$log_cac_z <- scale(data_expand$log_cac, center = TRUE, scale = TRUE)
data_expand$log_tac_z <- scale(data_expand$log_tac, center = TRUE, scale = TRUE)
data_expand$log_total_tac_z <- scale(data_expand$log_total_tac, center = TRUE, scale = TRUE)
data_expand$log_aortic_arch_ca_z <- scale(data_expand$log_aortic_arch_ca, center = TRUE, scale = TRUE)
data_expand$log_crp_mean_z <- scale(data_expand$log_crp_mean, center = TRUE, scale = TRUE)
data_expand$log_il1b_cytokines_mean_z <- scale(data_expand$log_il1b_cytokines_mean, center = TRUE, scale = TRUE)
data_expand$log_il6_cytokines_mean_z <- scale(data_expand$log_il6_cytokines_mean, center = TRUE, scale = TRUE)
data_expand$log_tnfa_cytokines_mean_z <- scale(data_expand$log_tnfa_cytokines_mean, center = TRUE, scale = TRUE)
data_expand$log_tareas_mean_z <- scale(data_expand$log_tareas_mean, center = TRUE, scale = TRUE)

# Bin CAC into categories
data_expand$cac_severity <- cut(
  data_expand$cac,
  breaks = c(-Inf, 0, 100, 400, Inf),
  labels = c("none", "mild", "moderate", "severe"),
  right = TRUE  # includes upper bound in interval, e.g. 100 is included in "mild"
)

# Combine moderate and severe CAC bins, because only one severe case
data_expand$cac_severity <- as.character(data_expand$cac_severity)  # convert to character for recoding
data_expand$cac_severity[data_expand$cac_severity %in% c("moderate", "severe")] <- "moderate/severe"
data_expand$cac_severity <- factor(data_expand$cac_severity, levels = c("none", "mild", "moderate/severe")) # Convert back to factor and re-level

# Bin TAC into four categories
data_expand$tac_severity <- cut(
  data_expand$tac,
  breaks = c(-Inf, 0, 100, 400, Inf),
  labels = c("none", "mild", "moderate", "severe"),
  right = TRUE  # includes upper bound in interval, e.g. 100 is included in "mild"
)

# Bin TAC into only three categories (merge moderate and severe)
data_expand$tac_severity_3bin <- cut(
  data_expand$tac,
  breaks = c(-Inf, 0, 100, Inf),
  labels = c("none", "mild", "moderate/severe"),
  right = TRUE  # includes upper bound in interval
)

# Bin AAC into four categories
data_expand$aac_severity <- cut(
  data_expand$aortic_arch_ca,
  breaks = c(-Inf, 0, 100, 400, Inf),
  labels = c("none", "mild", "moderate", "severe"),
  right = TRUE  # includes upper bound in interval, e.g. 100 is included in "mild"
)

# Bin AAC into only three categories (merge moderate and severe)
data_expand$aac_severity_3bin <- cut(
  data_expand$aortic_arch_ca,
  breaks = c(-Inf, 0, 100, Inf),
  labels = c("none", "mild", "moderate/severe"),
  right = TRUE  # includes upper bound in interval
)

# Create aortic calcification index = standardized mean of TAC and AAC
data_expand$aortic_calcification_index <- with(data_expand,
                                               ifelse(!is.na(log_tac_z) & !is.na(log_aortic_arch_ca_z),
                                                      (log_tac_z + log_aortic_arch_ca_z) / 2,
                                                      NA)  # assign NA if either value is missing
)

data_expand$global_calcification_burden <- with(data_expand,
                                                ifelse(!is.na(aortic_calcification_index) & !is.na(log_cac_z),
                                                       0.675 * aortic_calcification_index + 0.39 * log_cac_z,
                                                       NA)  # assign NA if either component is missing
)

# Create quartiles for number of tareas
data_expand <- data_expand %>%
  mutate(tareas_quartile = ntile(tareas_cuantos_mean, 4))

#prep data
all_dat <- data.frame (
  sex = data_expand$male,
  bmd_z = data_expand$bmd_z,
  age_z = data_expand$age_z,
  ffm_z = data_expand$ffm_z,
  thoracic_fat_z = data_expand$thoracic_fat_z,
  log_crp_mean_z = data_expand$log_crp_mean_z,
  log_il1b_cytokines_mean_z = data_expand$log_il1b_cytokines_mean_z,
  log_il6_cytokines_mean_z = data_expand$log_il6_cytokines_mean_z,
  log_tnfa_cytokines_mean_z = data_expand$log_tnfa_cytokines_mean_z,
  log_cac_z = data_expand$log_cac_z,
  log_tac_z = data_expand$log_tac_z,
  log_aortic_arch_ca_z = data_expand$log_aortic_arch_ca_z,
  log_cac = data_expand$log_cac,
  log_tac = data_expand$log_tac,
  log_aac = data_expand$log_aortic_arch_ca,
  cac = data_expand$cac,
  tac = data_expand$tac,
  aac = data_expand$aortic_arch_ca
)

dat <- make_dat(
  all_dat, 
  imput = TRUE
)

age_mean <- mean(data_expand$age)

# miss_inflam <- which(
#   is.na(all_dat$log_il1b_cytokines_mean_z) |
#     is.na(all_dat$log_il6_cytokines_mean_z)  |
#     is.na(all_dat$log_tnfa_cytokines_mean_z)
# )
# 
# df <- all_dat[complete.cases(all_dat), ]
# 
# dat <- list (
#   N = nrow(df),
#   S = 2,
#   sex = df$sex+1,
#   bmd_z = df$bmd_z,
#   age_z = df$age_z,
#   ffm_z = df$ffm_z,
#   thoracic_fat_z = df$thoracic_fat_z,
#   crp = df$log_crp_mean_z,
#   inflam_markers = as.matrix(data.frame( crp = df$log_crp_mean_z,
#                      il1b = df$log_il1b_cytokines_mean_z,
#                      il6 = df$log_il6_cytokines_mean_z,
#                      tfna = df$log_tnfa_cytokines_mean_z)),
#   #il1b = df$log_il1b_cytokines_mean_z,
#   #il6 = df$log_il6_cytokines_mean_z,
#   #tnfa = df$log_tnfa_cytokines_mean_z,
#   calc_markers = as.matrix(data.frame( cac = df$log_cac,
#                                        tac = df$log_tac,
#                                        aac = df$log_aac)),
#   any_calc = as.integer(df$log_cac > 0 |
#                        df$log_tac > 0 |
#                        df$log_aac > 0)
# )
