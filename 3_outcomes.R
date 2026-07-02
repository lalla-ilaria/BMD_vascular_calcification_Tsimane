
#extract samples
post_nosex <- lapply(fits_nosex, extract.samples)
post_c <- lapply(fits_c, extract.samples)
post_2 <- extract.samples(fits_2a)
post_2b <- extract.samples(fits_2b)
post_3 <- extract.samples(fits_3)
post_4pln <- lapply(fits_4_pzln, extract.samples)
post_5 <- lapply(fits_5, extract.samples)
post_6 <- lapply(fits_6, extract.samples)

#post_nosex_i <- lapply(fits_nosex_i, extract.samples)
#post_i <- lapply(fits_i, extract.samples)

#######################
# Calculate contrasts #
#######################

#defines calcification values to compare on the log scale
calc_0.1 <- log10(1.1)
#calc_100 <- log10(101)
calc_400 <- log10(401)

#Prep age for counterfactual construction
age_mean <- mean(data_expand$age)
age_SD <- sd(data_expand$age)
counterfactual_age_zscore <- (counterfactual_age - age_mean) / age_SD
counterfactual_age <- age_mean

gq <- list()



# No-sex model: BMD contrasts

gq <- gq_add(gq, name = "BMD_cont_calc0_0.1",
             interpretation = "mean and 89%PI of difference between BMD of average individual with no calcification and average individual with minimal calcification (0.1)",
             draws = counterfactual_BMD(post_nosex[[1]], calc = c(0,0,0)) - counterfactual_BMD(post_nosex[[1]], calc = c(calc_0.1,0,0)),
             section = "no_sex")

gq <- gq_add(gq, name = "BMD_slope_coeff_beta2",
             interpretation = "coefficient mean and PI for beta2 (i.e. the slope of calcification once there is any calcification) for a model with no sex differences, all calcification markers in the composite",
             draws = post_nosex[[1]]$beta2,
             section = "no_sex")

# by-marker contrasts, no-sex model, 0 to 400
# post_nosex[[2]] = CAC, post_nosex[[3]] = TAC, post_nosex[[4]] = AAC
nosex_markers <- list(cac = post_nosex[[2]], tac = post_nosex[[3]], aac = post_nosex[[4]])

for(marker_name in names(nosex_markers)){
  gq <- gq_add(gq, name = paste0("BMD_", marker_name, "_cont_calc0_400"),
               interpretation = paste0("mean and 89%PI of difference between BMD of average individual with no calcification and average individual with high risk calcification (400) in ", toupper(marker_name)),
               draws = counterfactual_BMD(nosex_markers[[marker_name]], calc = 0) - counterfactual_BMD(nosex_markers[[marker_name]], calc = calc_400),
               section = "no_sex")
}

for(marker_name in names(nosex_markers)){
  gq <- gq_add(gq, name = paste0("BMD_", marker_name, "_cont_calc0_0.1"),
               interpretation = paste0("mean and 89%PI of difference between BMD of average individual with no calcification and average individual with high risk calcification (400) in ", toupper(marker_name)),
               draws = counterfactual_BMD(nosex_markers[[marker_name]], calc = 0) - counterfactual_BMD(nosex_markers[[marker_name]], calc = calc_0.1),
               section = "no_sex")
}

# by-marker contrasts, no-sex model, 0.1 to 400
for(marker_name in names(nosex_markers)){
  gq <- gq_add(gq, name = paste0("BMD_", marker_name, "_cont_calc0.1_400"),
               interpretation = paste0("mean and 89%PI of difference between BMD of average individual with no calcification and average individual with high risk calcification (400) in ", toupper(marker_name)),
               draws = counterfactual_BMD(nosex_markers[[marker_name]], calc = calc_0.1) - counterfactual_BMD(nosex_markers[[marker_name]], calc = calc_400),
               section = "no_sex")
}

# Sex differences: BMD contrasts (no marker, post_c[[1]])

gq <- gq_add(gq, name = "BMD_men_cont_BMD_calc0_0.1",
             interpretation = "mean and 89%PI of difference between BMD of average individual with no calcification and average men with minimal calcification (0.1)",
             draws = counterfactual_BMD(post_c[[1]], sex = 2, calc = c(0,0,0)) - counterfactual_BMD(post_c[[1]], sex = 2, calc = c(calc_0.1,0,0)),
             section = "sex_stratified")

gq <- gq_add(gq, name = "BMD_women_cont_BMD_calc0_0.1",
             interpretation = "mean and 89%PI of difference between BMD of average woman with no calcification and average individual with minimal calcification (0.1)",
             draws = counterfactual_BMD(post_c[[1]], sex = 1, calc = c(0,0,0)) - counterfactual_BMD(post_c[[1]], sex = 1, calc = c(calc_0.1,0,0)),
             section = "sex_stratified")

gq <- gq_add(gq, name = "BMD_slope_coeff_beta2_men",
             interpretation = "coefficient mean and PI for beta2 (i.e. the slope of calcification once there is any calcification) for men, all calcification markers in the composite",
             draws = post_c[[1]]$beta2[,2],
             section = "sex_stratified")

gq <- gq_add(gq, name = "BMD_slope_coeff_beta2_women",
             interpretation = "coefficient mean and PI for beta2 (i.e. the slope of calcification once there is any calcification) for women, all calcification markers in the composite",
             draws = post_c[[1]]$beta2[,1],
             section = "sex_stratified")

# Sex-stratified, by-marker contrasts (post_c[[2]]=CAC, [[3]]=TAC, [[4]]=AAC)

sex_markers <- list(cac = post_c[[2]], tac = post_c[[3]], aac = post_c[[4]])
sex_labels  <- list(men = 2, women = 1) # sex index used in counterfactual_BMD()

# 0 to 0.1, by sex and marker
for(sex_name in names(sex_labels)){
  for(marker_name in names(sex_markers)){
    gq <- gq_add(gq, name = paste0("BMD_", sex_name, "_", marker_name, "_cont_calc0_0.1"),
                 interpretation = paste0("mean and 89%PI of difference between BMD of average individual with no calcification and average", sex_name, "with low calcification (0.1) in ", toupper(marker_name)),
                 draws = counterfactual_BMD(sex_markers[[marker_name]], sex = sex_labels[[sex_name]], calc = 0) - counterfactual_BMD(sex_markers[[marker_name]], sex = sex_labels[[sex_name]], calc = calc_0.1),
                 section = "sex_stratified")
  }
}

# 0.1 to 400, by sex and marker (only men + women had this in the original; aac/cac/tac all included)
for(sex_name in names(sex_labels)){
  for(marker_name in names(sex_markers)){
    gq <- gq_add(gq, name = paste0("BMD_", sex_name, "_", marker_name, "_cont_calc0.1_400"),
                 interpretation = paste0("mean and 89%PI of difference between BMD of average", sex_name, "with low calcification (0.1) and average individual with high risk calcification (400) in ", toupper(marker_name)),
                 draws = counterfactual_BMD(sex_markers[[marker_name]], sex = sex_labels[[sex_name]], calc = calc_0.1) - counterfactual_BMD(sex_markers[[marker_name]], sex = sex_labels[[sex_name]], calc = calc_400),
                 section = "sex_stratified")
  }
}

# 0 to 400, by sex and marker
for(sex_name in names(sex_labels)){
  for(marker_name in names(sex_markers)){
    gq <- gq_add(gq, name = paste0("BMD_", sex_name, "_", marker_name, "_cont_calc0_400"),
                 interpretation = paste0("mean and 89%PI of difference between BMD of average", sex_name, "with no calcification and average individual with high risk calcification (400) in ", toupper(marker_name)),
                 draws = counterfactual_BMD(sex_markers[[marker_name]], sex = sex_labels[[sex_name]], calc = 0) - counterfactual_BMD(sex_markers[[marker_name]], sex = sex_labels[[sex_name]], calc = calc_400),
                 section = "sex_stratified")
  }
}

doc <- read_docx()

for(i in seq_along(gq)){
  x <- gq[[i]]
  summ <- data.frame(Mean  = x$mean, Lower = x$PI[1], Upper = x$PI[2])
  
  doc <- body_add_par(doc, value = names(gq[i]) ,style = "heading 2")
  
  doc <- body_add_par(doc, value = paste("Section:", x$section), style = "Normal")
  
  doc <- body_add_par(doc, value = x$interpretation, style = "Normal")
  
  doc <- body_add_flextable(doc, flextable(summ))
  
  doc <- body_add_par(doc, "", style = "Normal")
}

print(doc, target = "outcomes/Generated_quantities.docx")
# Export

# str(gq, max.level = 2) # uncomment to inspect
# gq_export_json(gq, "gq.json") # then: node gq_to_word.js gq.json report.docx

######################################
# Figures
######################################

#########################
#FIGURE 2
#########################
#Figure 2 shows the association between calcification and BMD. All results are produced using linear regression models including the full set of covariates. In these models, calcification presence (i.e. 0 vs. any calcification) and amount are modeled separately (through parameters ????? and ?????). The resulting predictions are shown on each panel: the left side shows expected BMD for an individual with no calcification (mean and 89% interval of the posterior values calculated for an individual of average age, inflammation, thoracic fat and fat free mass). The right side shows hoe BMD is expected to decrease with increasing levels of calcification. Panels on the left column show results from models that do not explicitly model sex, while the models resulting in the right column plots allow sex to vary the slope for age and calcification. The top row shows results from a model where calcification is a composite of the three different vascular beds (calculated within the model using simplex weights, see SI section XXX and figure SI XXX for further explanation), while the following rows are relative to models where each vascular bed is analyzed separately (see SI section XXX for how these compare to the bed specific predictions for the composite model). 

# calcification variables
calc_vars <- c("cac", "tac", "aac")
letters_plot <- c("A", "E", "B", "F", "C", "G", "D", "H")

png("outcomes/Figure_2.png", width = 16, height = 24, units = "cm", res = 500, pointsize = 12)
#prepare plot layout
par(mfrow = c(4, 2),
    mar = c(3, 3, 2, 1),   # inner margins
    oma = c(0, 3, 0, 0),   # outer left margin for row labels
    mgp = c(1.8, 0.5, 0),
    tcl = -0.2
)
seqalong_c <- seq(0,8, length.out = 100)

covariates <- list(calc = c(0, 0, 0),
                   age = counterfactual_age_zscore,
                   ffm = 0,
                   thoracic_fat = 0,
                   inflam = c(0,0,0,0),
                   calc_composite = TRUE)

plot_counterfactual(
  post = post_nosex[["all"]],
  vary = "calc",
  vary_seq = seqalong_c,
  covariates = covariates,
  xlab = "Calcification composite",
  xlim = c(0,3.5),
  latent_points = FALSE,
  latent_var = "calc",
  dat = dat,
  col = allcol,
  main = "Sexes combined"
)
mtext(letters_plot[1], side = 3, line = 0, at = -0.9 )
mtext("composite", side = 2, line = 4,
      font = par("font.main"))

plot_counterfactual_sex(
  post = post_c[["all"]],
  vary = "calc",
  vary_seq = seqalong_c,
  covariates = covariates,
  xlab = "Calcification composite",
  xlim = c(0,3.5),
  latent_points = FALSE,
  latent_var = "calc",
  dat = dat,
  main = "Males vs. Females"
)
legend_ribbon(x0 = 2.2, y0 = 0.9,
  label = "Women", col = femalecol)
legend_ribbon(x0 = 2.2, y0 = 0.7,
  label = "Men", col = malecol)
mtext(letters_plot[2], side = 3, line = 0, at = -0.9 )

covariates$calc <- 0
for(v in calc_vars){
  idx <- match(v, calc_vars)
  
  plot_counterfactual(
    post = post_nosex[[v]],
    vary = "calc",
    vary_seq = seqalong_c,
    covariates = covariates,
    xlab = paste(v, "(log)"),
    xlim = c(0,3.5),
    col = allcol
  )
  mtext(v, side = 2, line = 4,
        font = par("font.main"))
  mtext(letters_plot[1+2*idx], side = 3, line = 0, at = -0.9 )
  #points(dat$calc_markers[,idx], dat$bmd_z, col = adjustcolor(cols[dat$sex], 0.3), pch = 16)
  
  plot_counterfactual_sex(
    post = post_c[[v]],
    vary = "calc",
    vary_seq = seqalong_c,
    covariates = covariates,
    xlab = paste(v, "(log)"),
    xlim = c(0,3.5)
  )
  mtext(letters_plot[2+2*idx], side = 3, line = 0, at = -0.9 )
  #points(dat$calc_markers[,idx], dat$bmd_z, col = adjustcolor(cols[dat$sex], 0.3), pch = 16)
}

dev.off()

##############################
# FIGURE 3
#############################

png("outcomes/Figure_3.png", width = 16, height = 6, units = "cm", res = 500, pointsize = 8)
par(mfrow = c(1, 2),
    mar = c(3, 3, 2, 1),   # inner margins
    oma = c(0, 0, 0, 0),   # outer left margin for row labels
    mgp = c(1.8, 0.5, 0),
    tcl = -0.2
)
covariates <- list(calc = c(0, 0, 0),
                   age = 0,
                   ffm = 0,
                   thoracic_fat = 0,
                   inflam = c(0,0,0,0),
                   inflam_composite = TRUE)

plot_counterfactual(
  post = post_nosex[[1]],
  vary = "inflam",
  vary_seq = seqalong,
  covariates = covariates,
  xlab = "Inflammation composite",
  latent_var = "inflam",
  col = directcol
)

plot_counterfactual(
  post = post_2,
  vary = "inflam",
  vary_seq = seqalong,
  covariates = covariates,
  xlab = "Inflammation composite",
  latent_var = "inflam",
  col = totalcol,
  add = T
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.8,
  label = "Direct (adj. for all covariates)",
  col = directcol
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.65,
  label = "Total (omitting calcification)",
  col = totalcol
)
mtext("A", side = 3, line = 0, at = -2 )

#coefficient posterior estimates
par(mar = c(3, 8, 2, 1))   # inner margins

post_inflammation_seq <- list(post_nosex[[1]], post_2)
plot(NULL, xlim = c(-0.5,0.5), ylim = c(0.5, 2.5), xlab = "Inflammation coefficient (??)", ylab = "", yaxt='n' )
abline(v = 0, lty = 2, col = "grey80")
for (i in 1:2){
  post <- post_inflammation_seq[[i]]
  points( mean(post$theta), i , pch = 16, cex = 2, col = c(directcol, totalcol)[i] )
  lines( PI(post$theta), rep(i, 2), lwd = 3, col = c(directcol, totalcol)[i])
}
axis(side = 2, at = 1:2, labels = c("Direct \n (adj. for all covariates)", "Total \n (omitting calcification)" ), las = 2)
mtext("B", side = 3, line = 0, at = -0.5 )
  
dev.off()

############################
# FIGURE 4
############################
png("outcomes/Figure_4.png", width = 16, height = 12, units = "cm", res = 500, pointsize = 11)
par(mfrow = c(2, 2),
    mar = c(3, 3, 2, 1),   # inner margins
    oma = c(0, 0, 0, 0),   # outer left margin for row labels
    mgp = c(1.8, 0.5, 0),
    tcl = -0.2
)

# Thoracic fat direct (all covariates) total (omitting calcification and inflammation)
covariates <- list(calc = c(0, 0, 0),
                   age = 0,
                   ffm = 0,
                   thoracic_fat = 0,
                   inflam = c(0,0,0,0))

plot_counterfactual(
  post = post_nosex[[1]],
  vary = "thoracic_fat",
  vary_seq = seqalong,
  covariates = covariates,
  xlab = "Thoracic fat  (z-scored)",
  col = directcol
)

plot_counterfactual(
  post = post_3,
  vary = "thoracic_fat",
  vary_seq = seqalong,
  covariates = covariates,
  col = totalcol,
  add = T
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.8,
  label = "Direct (adj. for all covariates)",
  col = directcol
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.55,
  label = "Total (omitting calcification \n and inflammation)",
  col = totalcol
)
mtext("A", side = 3, line = 0, at = -2 )


#coefficient posterior estimates
par(mar = c(3, 8, 2, 1))   # inner margins

post_inflammation_seq <- list(post_nosex[[1]], post_2, post_2b, post_3)
plot(NULL, xlim = c(-0.5,0.5), ylim = c(0.5, 4.5), xlab = "Thoracic fat coefficient (??)", ylab = "", yaxt='n' )
abline(v = 0, lty = 2, col = "grey80")
for (i in 1:4){
  post <- post_inflammation_seq[[i]]
  points( mean(post$eta), i , pch = 16, cex = 2, col = c(directcol,partialcol1, partialcol2, totalcol)[i] )
  lines( PI(post$eta), rep(i, 2), lwd = 3, col = c(directcol,partialcol1, partialcol2, totalcol)[i])
}
axis(side = 2, at = 1:4, labels = c("Direct \n (adj. for all covariates)", "Through calcification", "Through inflammation", "Total \n (omitting calcification \n and inflammation)" ), las = 2)
mtext("D", side = 3, line = 0, at = -0.5 )

par(mar = c(3, 3, 2, 1))   # inner margins

# Thoracic fat direct (all covariates) through calc (omitting calcification)
plot_counterfactual(
  post = post_nosex[[1]],
  vary = "thoracic_fat",
  vary_seq = seqalong,
  covariates = covariates,
  xlab = "Thoracic fat  (z-scored)",
  col = directcol
)

plot_counterfactual(
  post = post_2,
  vary = "thoracic_fat",
  vary_seq = seqalong,
  covariates = covariates,
  col = partialcol1,
  add = T
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.8,
  label = "Direct (adj. for all covariates)",
  col = directcol
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.65,
  label = "Through calcification",
  col = partialcol1
)
mtext("B", side = 3, line = 0, at = -2 )

# Thoracic fat direct (all covariates) through infl (omitting inflammation)
plot_counterfactual(
  post = post_nosex[[1]],
  vary = "thoracic_fat",
  vary_seq = seqalong,
  covariates = covariates,
  xlab = "Thoracic fat  (z-scored)",
  col = directcol
)

plot_counterfactual(
  post = post_2b,
  vary = "thoracic_fat",
  vary_seq = seqalong,
  covariates = covariates,
  col = partialcol2,
  add = T
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.8,
  label = "Direct (adj. for all covariates)",
  col = directcol
)

legend_ribbon(
  x0 = -1.5,
  y0 = 0.65,
  label = "Through inflammation",
  col = partialcol2
)
mtext("C", side = 3, line = 0, at = -2 )
dev.off()




####################
# Figure 5
####################
png("outcomes/Figure_5.png", width = 16, height = 12, units = "cm", res = 500, pointsize = 8)
par(mfrow = c(2,2),
    mar = c(4,4,2,1)
)
inflam_vars <- c("crp","il1b","il6","tnfa")
inflam_varsC <- c("CRP","IL-1??","IL-6","TNF-??")
covariates <- list(age = counterfactual_age_zscore,
                   ffm = 0,
                   thoracic_fat = 0)
for(v in inflam_vars){
  idx <- match(v, inflam_vars)
  #thoracic_fat
  plot_counterfactual(
    post = post_6[[v]],
    outcome = "infl",
    vary = "thoracic_fat",
    xlab = "Thoracic fat  (z-scored)",
    ylab = paste("Expected", inflam_varsC[idx], "(z-scored)"),
    xlim = c(-2, 2),
    ylim = ylim, 
    main = inflam_varsC[idx]
  )
  mtext(LETTERS[idx], side = 3, line = 0, at = -2 )
  

}
dev.off()
#########################################################################################################
# # Inflammation variables
# inflam_vars <- c("crp","il1b","il6","tfna")
# seqalong <- seq(-4,4, length.out = 100)
# #letters_plot <- c("A", "F", "B", "G", "C", "H", "D", "I", "E", "K")
# 
# png("outcomes/Figure_5.png", width = 8, height = 6, units = "cm", res = 500, pointsize = 12)
# #prepare plot layout
# par(mfrow = c(1, 1),
#     mar = c(3, 3, 2, 1),   # inner margins
#     oma = c(0, 0, 0, 0),   # outer left margin for row labels
#     mgp = c(1.8, 0.5, 0),
#     tcl = -0.2
# )
# 
# covariates <- list(calc = c(0, 0, 0),
#                    age = 0,
#                    ffm = 0,
#                    thoracic_fat = 0,
#                    inflam = c(0,0,0,0),
#                    inflam_composite = TRUE)
# 
# plot_counterfactual(
#   post = post_nosex[["all"]],
#   vary = "inflam",
#   vary_seq = seqalong_c,
#   covariates = covariates,
#   xlab = "Inflammation composite",
#   xlim = c(0,3.5),
#   latent_points = FALSE,
#   latent_var = "infl",
#   dat = dat,
#   col = allcol,
#   main = "Sexes combined"
# )
# mtext("A", side = 3, line = 0, at = -4 )
# mtext("composite", side = 2, line = 4,
#       font = par("font.main"))


# covariates$calc <- 0
# for(v in inflam_vars){
#   idx <- match(v, inflam_vars)
#   
#   plot_counterfactual(
#     post = post_nosex_i[[v]],
#     vary = "inflam",
#     vary_seq = seqalong,
#     covariates = covariates,
#     xlab = paste(v, "(z-scored)"),
#     col = allcol,
#     main = v
#   )
#   
#   mtext(LETTERS[idx], side = 3, line = 0, at = -2 )
  #points(dat$infl_markers[,idx], dat$bmd_z, col = adjustcolor(cols[dat$sex], 0.3), pch = 16)
  # 
  # plot_counterfactual_sex(
  #   post = post_i[[v]],
  #   vary = "inflam",
  #   vary_seq = seqalong_c,
  #   covariates = covariates,
  #   xlab = paste(v, "(log)"),
  #   xlim = c(0,3.5)
  # )
  # mtext(letters_plot[2+2*idx], side = 3, line = 0, at = -4 )
  #points(dat$infl_markers[,idx], dat$bmd_z, col = adjustcolor(cols[dat$sex], 0.3), pch = 16)
#}

#dev.off()

# ######################################################################################################################
# covariates <- list(calc = c(0, 0, 0),
#                    age = 0,
#                    ffm = 0,
#                    thoracic_fat = 0,
#                    inflam = c(0,0,0,0),
#                    calc_composite = TRUE)
# 
# post_nosex <- extract.samples(fits_nosex)
# d_nosex <- counterfactual_line(
#   post = post_nosex,
#   sex = 1,
#   covariates = covariates,
#   vary = "calc",
#   vary_seq = seq(0,4, length.out = 100)
# )
# 
# post_c <- lapply(fits_c, extract.samples)
# post_imput <- extract.samples(fits_imput)
# 
# png("outcomes/composites.png", width = 24, height = 14, units = "cm", res = 500)
# par(mfrow = c(2,1), mar = mar, mgp = c(2, 1.2, 0))
# covariates <- list(calc = c(0, 0, 0),
#                    age = 0,
#                    ffm = 0,
#                    thoracic_fat = 0,
#                    inflam = c(0,0,0,0),
#                    calc_composite = TRUE)
# 
# plot_counterfactual_sex(
#   post = post_c[["all"]],
#   vary = "calc",
#   vary_seq = seq(0,4, length.out = 100),
#   covariates = covariates,
#   xlab = "Calcification composite",
#   xlim = c(0,3.5),
#   latent_points = FALSE,
#   latent_var = "calc",
#   dat = dat,
#   main = paste("With complete cases only. Sample size =", dat$N)
# )
# 
# plot_counterfactual_sex(
#   post = post_imput,
#   vary = "calc",
#   vary_seq = seq(0,4, length.out = 100),
#   covariates = covariates,
#   xlab = "Calcification composite",
#   xlim = c(0,3.5),
#   latent_points = FALSE,
#   latent_var = "calc",
#   dat = dat_imput,
#   main = paste("With imput of missing data. Sample size =", dat_imput$N)
# )
# 
# plot_counterfactual(
#   post = post_nosex,
#   vary = "calc",
#   vary_seq = seq(0,4, length.out = 100),
#   covariates = covariates,
#   xlab = "Calcification composite",
#   xlim = c(0,3.5),
#   latent_points = FALSE,
#   latent_var = "calc",
#   dat = dat_nosex,
#   col = allcol,
#   main = paste("With complete cases, nosex. Sample size =", dat_nosex$N)
# )
# 
# covariates <- list(calc = c(0, 0, 0),
#                    age = 0,
#                    ffm = 0,
#                    thoracic_fat = 0,
#                    inflam = c(0,0,0,0),
#                    inflam_composite = TRUE)
# 
# 
# plot_counterfactual_sex(
#   post = post_c[["all"]],
#   vary = "inflam",
#   vary_seq = seqalong,
#   covariates = covariates,
#   xlab = "Inflammation composite",
#   latent_points = FALSE,
#   latent_var = "inflam",
#   dat = dat
# )
# 
# plot_counterfactual_sex(
#   post = post_imput,
#   vary = "inflam",
#   vary_seq = seqalong,
#   covariates = covariates,
#   xlab = "Inflammation composite",
#   latent_points = FALSE,
#   latent_var = "inflam",
#   dat = dat_imput
# )
# 
# plot_counterfactual(
#   post = post_nosex,
#   vary = "inflam",
#   vary_seq = seqalong,
#   covariates = covariates,
#   xlab = "Inflammation composite",
#   latent_points = FALSE,
#   latent_var = "inflam",
#   dat = dat_nosex,
#   col = allcol
# )
# 
# dev.off()
# 
