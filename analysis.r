## Analyses for Madin et al. GBR satellite imagery survey grazing halos paper

# SECTION 1: halo presence (as response variable)

# question: are halos more (or less) PREVALENT in reserves (or as reserves mature), areas that are more  
# heavily fished, and/or have diff't Chl a & temperature regimes?

# method: halo presence/absence as a function of all relevant predictor variables

## load data

datx <- read.csv("data/dat_public_halo_presence.csv", as.is=TRUE)

## explore

# check for correlations among variables (OK to incl both in model if R2 is < ~0.6)

pairs(datx[c("reserve_age", "temp_sst","chl_log","log_distance_mainland_km","abs_latitude")])

cor.test(datx$temp_sst, datx$abs_latitude) # Tight correl (0.94); keep temp_sst b/c more meaningful
cor.test(datx$log_distance_mainland_km, datx$temp_sst) # weak correl 
cor.test(datx$log_distance_mainland_km, datx$chl_log) # weak (0.10) correl. 

## analysis

# run single-variable models

# note: halos_present_bin is binary (0 vs 1) representation of halos_present
# compare all AIC values resuling from the various predictor variables below to AIC 
# from random model @ top (269.57) 

mod_halopresence_null <- glm(halos_present_bin ~ 1, family = binomial, data=datx) # 
summary(mod_halopresence_null)
# AIC: 269.57

mod_halopresence_reserve_age <- glm(halos_present_bin ~ reserve_age, family = binomial, data=datx) # 
summary(mod_halopresence_reserve_age)
# *, AIC: 266.53

mod_halopresence_zone <- glm(halos_present_bin ~ zone_broad_adj, family = binomial, data=datx) 
summary(mod_halopresence_zone)          # zone_broad_adj only counts =<8 year old reserves as protected; Babcock 2010
# *, AIC: 266.16

mod_halopresence_latitude <- glm(halos_present_bin ~ abs_latitude, family = binomial, data=datx) 
summary(mod_halopresence_latitude)
# ., AIC: 267.68

mod_halopresence_sst <- glm(halos_present_bin ~ temp_sst, family = binomial, data=datx) 
summary(mod_halopresence_sst)
hist(resid(mod_halopresence_sst))
# ., AIC: 268.43

mod_halopresence_chl <- glm(halos_present_bin ~ chl_log, family = binomial, data=datx) 
summary(mod_halopresence_chl)
hist(resid(mod_halopresence_chl))
# NS, AIC: 271.18

mod_halopresence_season <- glm(halos_present_bin ~ image_season, family = binomial, data=datx) 
summary(mod_halopresence_season)
anova(mod_halopresence_season, test="Chisq")
# NS, AIC: 239.11

mod_halopresence_mainlanddist <- glm(halos_present_bin ~ log_distance_mainland_km, family = binomial, data=datx) 
plot(jitter(halos_present_bin, 0.1) ~ log_distance_mainland_km, data=datx)
summary(mod_halopresence_mainlanddist)
# NS, AIC: 254.7

# drop variables

# based on the above correlations and single variable tests, we remove:
# zone_broad_adj b/c had to choose btwn that and reserve_age, and latter is more nuanced metric of same thing
# image_season because of data loss (N goes from 1300 to 377)
# latitude b/c tight correlation with temperature, and temp is more biologically meaningful
# log_distance_mainland_km b/c strong correl in halo size analyses w/temp_sst (even though not
# very strong correl in this dataset), so assuming both should not be included here (?)

# run models using all remaining variables

datx_reservesonly <- datx[datx$zone_broad=="protected",]
mod_halopresence_allvars_reservesonly <- glm(halos_present_bin ~ reserve_age + temp_sst + chl_log, family = binomial, data=datx_reservesonly) 
summary(mod_halopresence_allvars_reservesonly)
drop1(mod_halopresence_allvars_reservesonly, test="Chisq") 
step(mod_halopresence_allvars_reservesonly, test="Chisq") 

mod_halopresence_reservesonly_drop1 <- glm(halos_present_bin ~ reserve_age, family = binomial, data=datx_reservesonly)                                            # FINAL MODEL
summary(mod_halopresence_reservesonly_drop1)

# look at variance explained by remaining variables

library(hier.part)

datx_reservesonly_temp <- na.omit(datx_reservesonly[c("halos_present_bin", "reserve_age", "temp_sst", "chl_log")])                               # remove NAs, as required by the hier.part function
vars <- datx_reservesonly_temp[c("reserve_age", "temp_sst", "chl_log")]
hier.part(datx_reservesonly_temp$halos_present_bin, vars)                       

## plots

# manuscript fig. 2 (fig 1 is conceptual & made elsewhere)

## panel A: proportion of reefs with halos in fished vs protected reserves 
# (i.e., how many whole [not patch] reefs, where the habitat exists and we have clear imagery, have halos)

par(mfrow=c(1,2), mar=c(5, 4, 1, 3), oma=c(3, 1, 1, 1))

mod_halopresence_zone_specific <- glm(halos_present_bin ~ zone_specific, family = binomial, data=datx) 
summary(mod_halopresence_zone_specific)          
prediction_zone_specific <- predict(mod_halopresence_zone_specific, list(zone_specific = c("fished", "young reserve", "mature reserve")), se.fit=TRUE, type="response")
# type = "response" puts output back into original format (ie, proportion in this case) 
# se.fit=TRUE adds SE to model's output (default is SE not included)
plot(1:3, prediction_zone_specific$fit, ylim=c(0, 1), xlim=c(0.5, 3.5), axes=FALSE, ylab="Probability of halo occurrence", xlab="", pch = 19)
# 1:2 tells plot that x values are 1 and 2 (ie, x1 = 1 and x2 = 2)
# prediction_zone$fit has 2 values, so plot uses these as y1 & y2
# so we're giving plot 2 x vals & 2 y vals, which it then plots as normal
# pch says use filled circles as points
arrows(1:3, prediction_zone_specific$fit + prediction_zone_specific$se.fit, 1:3, prediction_zone_specific$fit - prediction_zone_specific$se.fit, code=3, angle=90, length=0.1)
# arrows are just an easy way to make error bars; takes same x/y input as plot, but repeats 
# to do the y values + prediction_zone$se.fit & then y values - prediction_zone$se.fit
# code=3 draws a head at both ends of the arrow 
# angle=90 flattens arrows so they're just lines, not points
# length=0.1 shortens arrow head (line) width from wider default 
axis(1, at=c(1, 2, 3), labels=c("Fished", "Young
                                reserve", "Mature
                                reserve"), las=2)
axis(2, las=2)  # las=2 orients y-axis ticks horizontally

## panel B: halo prevalence vs reserve age

mod_halopresence_reservesonly_drop1 <- glm(halos_present_bin ~ reserve_age, family = binomial, data=datx_reservesonly)           # just copying from above to have handy

library(fmsb)
NagelkerkeR2(mod_halopresence_reservesonly_drop1) # Calc R2 (goodness of fit) for glm model (Nagelkerke's R squared = model's power of explanation)
reserve_age_vec <- seq(min(datx_reservesonly$reserve_age, na.rm=TRUE), max(datx_reservesonly$reserve_age, na.rm=TRUE), 1)
prediction_reserve_age <- predict(mod_halopresence_reservesonly_drop1, list(reserve_age=reserve_age_vec, sqrt_catch_piscivores_nosharks=rep(0, length(reserve_age_vec))), type="response", se.fit=TRUE)

plot((reserve_age_vec/38)*11, prediction_reserve_age$fit, type="l", xlim=c(0, 12), ylim=c(0,1), axes=FALSE, xlab="", ylab="Probability of halo occurence", col="grey40")
title(xlab="Reserve age 
      (years)", line=3.5)           # this moves xlab down (before it was overlapping w/numbers)
polygon(c((reserve_age_vec/38)*11, rev((reserve_age_vec/38)*11)), c(prediction_reserve_age$fit+prediction_reserve_age$se.fit, rev(prediction_reserve_age$fit-prediction_reserve_age$se.fit)), col=rgb(0,0,0,0.3),border=NA)
axis(2, las=2)
axis(1, at=(seq(0, 40, 10)/38)*11, labels=seq(0, 40, 10), las=2)
points(jitter(((datx$reserve_age/38)*11), 10), jitter(((datx$halos_present_bin)), 0.1), col="black")


# SECTION 2: halo size (as response variable)

# question: are halos larger (or smaller) in reserve vs. non-reserve reefs (AND/OR as a function of latitude, patch area, etc.)?

# method: halo width as a function of all relevant predictor variables

## load data

datx <- read.csv("data/dat_public_halo_size.csv", as.is=TRUE)

## explore

# check for correlations among variables (OK to incl both in model if R2 is < ~0.6)

pairs(datx[c("log10_patch_area_m2","reserve_age", "temp_sst","chl_log","log_distance_mainland_km","abs_latitude")])

cor.test(datx$temp_sst, datx$abs_latitude)              # keep temp_sst b/c more biologically meaningful
cor.test(datx$chl_log, datx$log_distance_mainland_km)   # (0.26) could keep both 
cor.test(datx$temp_sst, datx$log_distance_mainland_km)  # drop log_distance_mainland_km; 
# (-0.73) sst more transferrable/bilolog'lly meaningful
cor.test(datx$temp_sst, datx$chl_log)                   # 0.43 - probably okay, but still some correl
cor.test(datx$abs_latitude, datx$log_distance_mainland_km) # (0.77) have already dropped both so no problem

cor.test(datx$log10_patch_area_m2, datx$log_distance_mainland_km) # no real correl.
cor.test(datx$log10_patch_area_m2, datx$abs_latitude) # (-0.37) some correl, but not huge
cor.test(datx$log10_patch_area_m2, datx$temp_sst) # (0.38) some correl, but not huge; same process as latitude
cor.test(datx$log10_patch_area_m2, datx$chl_log) # (-0.36) some correl, but not huge
cor.test(datx$log10_patch_area_m2, datx$reserve_age) # no real correl.

boxplot(datx$log_distance_mainland_km, as.factor(datx$zone_broad)) 
dev.off()

# do geometric tests to see what kind of process is governing halo size 

# (eg, is it a function of patch reef area, perimeter, rugosity...this should tell us something about 
# what biological process is governing halo size.)
# ANS: halo area is more a function of patch reef area than patch reef perimeter, based 
# on which one had lowest AIC  

mod_random <- lm(log10_halo_width ~ 1, data=datx)  # random model to base comparisons against
summary(mod_random)
AIC(mod_random)

mod_patch_perim <- lm(log10_halo_width ~ log10_patch_perimeter_m, data=datx)   # perimeter model
summary(mod_patch_perim)
AIC(mod_patch_perim)

mod_patch_area <- lm(log10_halo_width ~ log10_patch_area_m2, data=datx)  # patch area is best geometric model
summary(mod_patch_area)
AIC(mod_patch_area)

plot(log10_halo_width ~ log10_patch_area_m2, data=datx, col="grey")
ss <- seq(min(datx$log10_patch_area_m2, na.rm=TRUE), max(datx$log10_patch_area_m2, na.rm=TRUE), 0.1)
predict_mod_patch_area <- predict(mod_patch_area, list(log10_patch_area_m2=ss), interval="confidence")
lines(ss, predict_mod_patch_area[,"fit"])
lines(ss, predict_mod_patch_area[,"upr"], lty=2)
lines(ss, predict_mod_patch_area[,"lwr"], lty=2)

## analysis

# run most basic model
# compare all AIC values resulting from the various predictor variables below to AIC from random model @ top
mod_halosize_null <- lm(log10_halo_width ~ 1, data=datx) # needs family = binomial?
summary(mod_halosize_null)
AIC(mod_halosize_null)
# AIC: 58.85466

mod_halosize_reserve_age <- lm(log10_halo_width ~ reserve_age, data=datx) 
summary(mod_halosize_reserve_age)
AIC(mod_halosize_reserve_age)
# *; AIC: 56.38727 

mod_halosize_zone_broad <- lm(log10_halo_width ~ zone_broad, data=datx) 
summary(mod_halosize_zone_broad)            
AIC(mod_halosize_zone_broad)
# NS; AIC: 60.15179 

mod_halosize_zone_broad_adj <- lm(log10_halo_width ~ zone_broad_adj, data=datx) 
summary(mod_halosize_zone_broad_adj)        # zone_broad_adj only counts <=8 year old reserves as protected
AIC(mod_halosize_zone_broad_adj)
# **; AIC: 57.80198 

mod_halosize_sst <- lm(log10_halo_width ~ temp_sst, data=datx) 
summary(mod_halosize_sst)
AIC(mod_halosize_sst)
# ***; AIC: -39.84485  

mod_halosize_chl <- lm(log10_halo_width ~ chl_log, data=datx) 
summary(mod_halosize_chl)
AIC(mod_halosize_chl)
# ***; AIC: -25.72507

mod_halosize_mainlanddist <- lm(log10_halo_width ~ log_distance_mainland_km, data=datx) 
summary(mod_halosize_mainlanddist)
AIC(mod_halosize_mainlanddist)
# ***; AIC: 16.25269

# drop variables

# based on the above correlations and single variable tests, we remove:
# zone_broad_adj b/c had to choose btwn that and reserve_age, and latter is more nuanced metric of same thing
# image_season because of data loss 
# latitude b/c tight correlation with temperature, and temp is more biologically meaningful
# log_distance_mainland_km b/c strong correl in halo size analyses w/temp_sst (even though not
# very strong correl in this dataset), so assuming both should not be included here (?)

# run models using all remaining variables

# regular linear model
datx_reservesonly <- datx[datx$zone_broad=="protected",]
mod_halowidth_allvars_reservesonly <- lm(log10_halo_width ~ log10_patch_area_m2 + reserve_age + temp_sst + chl_log, data=datx_reservesonly) 

# mod_halowidth_allvars_reservesonly <- lm(log10_halo_width ~ reserve_age + temp_sst + chl_log + log_distance_mainland_km, data=datx_reservesonly)                         # adding in log_distance_mainland_km
summary(mod_halowidth_allvars_reservesonly)
drop1(mod_halowidth_allvars_reservesonly, test="F") 
step(mod_halowidth_allvars_reservesonly, test="F") 

# regular linear model, but with interaction term 
mod_halowidth_reservesonly_drop1_intrctn <- lm(log10_halo_width ~ log10_patch_area_m2 * temp_sst, data=datx_reservesonly)                               # check to see if interactions signif  
summary(mod_halowidth_reservesonly_drop1_intrctn) # says there is a signif intrctn, so investigate with coplot
with(datx_reservesonly, coplot(log10_halo_width ~ temp_sst | log10_patch_area_m2, panel=panel.smooth))
# shows that the intrctn just driven by a few low temp points,
# so prob not substantial enough to include in final model 
# 'with' just means don't have to attach data with $ each time
mod_halowidth_reservesonly_drop1 <- lm(log10_halo_width ~ log10_patch_area_m2 + temp_sst, data=datx_reservesonly)                                           # final regular linear model, but FINAL model is below
summary(mod_halowidth_reservesonly_drop1)

# mixed-effects model (to incorp reef_id as random factor) - FINAL, CORRECT MODEL TO USE
    # since the spatially-structured nature of the data means that we could possibly have strong  
    # spatial pseudo-replication that we have to deal with 
    # ANS: this says that only patch area is significant...which model to use?
library(nlme)
mod_halowidth_reservesonly_lme <- lme(log10_halo_width ~ log10_patch_area_m2 + reserve_age + temp_sst + chl_log, data=datx_reservesonly, random = ~ 1 + log10_patch_area_m2 | reef_id, na.action = na.omit, method="ML") 
# Drew's suggested model: patch area is both fixed and
# random factor; latter allows intercept to vary, not just slope, as fnctn of reef_id

summary(mod_halowidth_reservesonly_lme)
names(summary(mod_halowidth_reservesonly_lme))    # lists elements that can be extracted from summary
drop1(mod_halowidth_reservesonly_lme, test="Chisq") # says only reef area is signif.

# also testing for interactions using final model
#mod_halowidth_reservesonly_lme_intcn <- lme(log10_halo_width ~ log10_patch_area_m2 * reserve_age * temp_sst * chl_log, data=datx_reservesonly, random = ~ 1 + log10_patch_area_m2 | reef_id, na.action = na.omit, method="ML") 
#summary(mod_halowidth_reservesonly_lme_intcn)

# look at variance explained by remaining variables

library(hier.part)

# note: hier.part function doesn't work for GLMM, so what Luiz and co.
# did in PNAS paper was presented the hier.part for the GLM (i.e., the regular model without random
# factors), not GLMM, model, and stated so in the text.

datx_reservesonly_temp <- na.omit(datx_reservesonly[c("log10_halo_width", "log10_patch_area_m2", "reserve_age", "temp_sst", "chl_log")])                        # remove NAs, as required by the hier.part function
vars <- datx_reservesonly_temp[c("log10_patch_area_m2", "reserve_age", "temp_sst", "chl_log")]
hier.part(datx_reservesonly_temp$log10_halo_width, vars)                       


## plots

# manuscript fig. 3 

datx$log10_patch_area_m2 <- log10(datx$patch_area_m2)
mod_halowidth_patch_area <- lm(log10_halo_width ~ log10_patch_area_m2, datx)

plot(log10(datx$patch_area_m2), datx$log10_halo_width, col="lightgrey", axes=FALSE, ylim=c(0, 2), xlab="", ylab="")  
# leave axes off so can more easily manipulate below
title(xlab="log10(Patch reef area), m2")      # add x-axis title
title(ylab="log10(Halo width), m")      # add y-axis title
axis(1, at=c(0, 1, 2, 3, 4), labels=10^c(0, 1, 2, 3, 4), las=2) # axis 1 is x; "at" is where to put axis labels;                                               # labels puts into non-log for readability; las is orientation
axis(2, at=c(0, 1, 2), labels=10^c(0, 1, 2), las=2) # axis 2 is y
# add model fit (prediction) & confidence intervals:
newx <- seq(min(datx$log10_patch_area_m2, na.rm = TRUE), max(datx$log10_patch_area_m2, na.rm = TRUE), length.out=50)                                          # makes 50-long seq of x-vals for predicted vals to plot against
# of same range as actual x-vals (patch reef areas)
predictedy <- predict(mod_halowidth_patch_area, newdata = data.frame(log10_patch_area_m2=newx), 
                      interval = 'confidence')     # make predicted y vals based on substituting newx for 
# log10_patch_area_m2; return conf intervals (can do SE, etc.)
polygon(c(rev(newx), newx), c(rev(predictedy[ ,3]), predictedy[ ,2]), col = rgb(0,0,0,0.2), border = NA)
# plot poylgon for x1, x2, 1, y2
# x1 = (rev(newx); x2 = newx; 
# y1 = rev(predictedy[ ,3] = upr CI; y2 = predictedy[ ,2] = lwr CI
lines(newx, predictedy[ ,1]) 
#text(3.5, 0.5, paste0("Slope=", round(mod_halowidth_patch_area$coefficients[2], 2), " (±", round(mod_halowidth_patch_area$coefficients[2] - confint(mod_halowidth_patch_area)[2,1], 2), " CIs)"))  # add text for slope/CIs

# add null model slopes for comparison:
origin_y <- predict(mod_halowidth_patch_area, list(log10_patch_area_m2=newx[length(newx)/2])) # sets origin as mod fit mdpt;
# can't start at 0 bc on log scale can't have 0, so would have
# to choose arbitrary origin, which is confusing 
origin_x <- newx[length(newx)/2]              # sets origin as midpoint of model fit
segments(origin_x, origin_y, newx[length(newx)], origin_y, lty=2) # adds flat line to represent slope=0, 
                                              # which is what we got from null model based on
                                              # purely perimeter-based geometric expectation 
segments(origin_x, origin_y, newx[length(newx)], newx[length(newx)] * 0.5, lty=2) # adds 0.5 slope line, 
                                              # which is what we got from null model based on  
                                              # purely area-based geometric expectation 


