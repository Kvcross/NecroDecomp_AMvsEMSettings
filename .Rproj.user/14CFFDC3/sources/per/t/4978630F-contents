# title: Summer 2017 Fungal Necromass Collaboration Mass Remaining Stats
# author: KVB
 #date: "9/21/2018"

# housekeeping 
rm(list=ls())

# load required packages
install.packages("pacman") # a package which allows you to load a bunch of packages at once
pacman::p_load(tidyverse,dplyr,broom,lmerTest,minpack.lm,Hmisc,car,gtools,lme4,nlme,minpack.lm,visreg)

# read in data file
fundecomp <- read.delim("data/CompiledNecroDecomp_summer2017.csv", sep = ",")

# specifying factors for data file
fundecomp  <-  within(fundecomp , {
            site = factor(site, ordered=FALSE)
            veg_type = factor(veg_type, ordered=FALSE)
            exp = factor(exp, ordered=FALSE)
            location = factor(location, ordered=FALSE)
            plot_Id = factor(plot_Id, ordered=FALSE)
            unique_plotID = factor(unique_plotID, ordered=FALSE)
            treatment = factor(treatment, ordered=FALSE)
            myc_assoc  = factor(myc_assoc, ordered=FALSE)
            isolate = factor(isolate, ordered=FALSE)
})

# removing NAs
fundecomp <- fundecomp%>% 
            filter(!is.na(mass_remaining))

sum(is.na(fundecomp$mt.mo))

# subsetting

fundecomp_myc <- fundecomp  %>%
  filter(myc_assoc != "mixed") 

fundecomp_melb.mort <- fundecomp  %>%
  filter(isolate !="Ceno" & isolate !="Mel_White") 
  
# averages 
mass.avg <- fundecomp %>%
            select(veg_type,isolate, mass_remaining) %>%
            group_by(veg_type,isolate) %>%
            summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(mass.avg , "output/avgs_massdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# data exploration
str(fundecomp)
summary(fundecomp)
# count the data to confirm that the design is unbalanced 
ftable(xtabs(~ isolate + incub_period.days + veg_type, fundecomp_melb.mort))

# specifying if factors are nested
with(fundecomp, isNested(veg_type,site))
with(fundecomp, isNested(incub_period.days,exp)) # isolate is not nested within veg type,but according to the cross tabulation isolate and incubation period are in partially crossed configurations
with(fundecomp, isNested(isolate,exp))

# looking at the distribution of mass remaining (mt.mo)
hist(fundecomp$mt.mo)

# looks a little right skewed- trying a log transformation
mt.mo1 <-  fundecomp$mt.mo + 1 # adding 1 to mt.mo to prevent negative values during log transformation
fundecomp$logmt.mo <- log(mt.mo1)
range(fundecomp$logmt.mo)
hist(fundecomp$logmt.mo)

#still looking a little skewed- now trying a logit transformation for proportion data
fundecomp$lmt.mo = (qlogis(fundecomp$mt.mo+0.001))  #lmt.mo is logit transformation of the percent disturbed plus 0.001 to prevent a zero proportion
range(fundecomp$lmt.mo)
hist(fundecomp$lmt.mo)
# looks more normal

# performing a logit transformation for datasheet which exlcudes Ceno and MB_white 
fundecomp_melb.mort$lmt.mo = (qlogis(fundecomp_melb.mort$mt.mo+0.001)) 

# doing a logistic transformation and then modelling linearly 
ggplot(fundecomp_melb.mort, aes(x = incub_period.days, y = lmt.mo, color=mel_status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Logit(Mass Remaining") +
  xlab("Days of field Incubation)") +
  facet_grid(. ~ veg_type)

ggplot(fundecomp_melb.mort, aes(x = incub_period.days, y = lmt.mo, color=mel_status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Logit(Mass Remaining") +
  xlab("Days of field Incubation)") +
  facet_wrap( veg_type) 

# looking at relationship with pH
ggplot(fundecomp_melb.mort, aes(x = pH, y = lmt.mo, color=isolate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Logit(Mass Remaining") +
  xlab("pH)") 

# looking for collinearity
vars <- c("soil_moisture", "pH", "incub_period.days","mt.mo")
pairs(fundecomp[,vars])
plot(fundecomp$soil_moisture ~ fundecomp$pH) 
cor.test(fundecomp$soil_moisture, fundecomp$pH) #low 0.3 collinearity between soil_moisture and pH


# using a linear mixed effect model to look at the effects of incubation period, myc association, veg_type and fungal isolate (fixed) and  
# I will fit the most complex model and then eliminate fixed effects to determine if the eliminated effect is significant
# the experimental design is partially crossed as we don't have data for each isolate in each vegetation type for each sampling date
# to account for the partially crossed design I will include a random effects which groups incubation periods within isolates/ experiment
# full model
mrlmer_full <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days +  myc_assoc + soil_moisture:pH + myc_assoc:pH +
                       (1|incub_period.days:(isolate:exp)) ,data = fundecomp,  REML=FALSE)
summary(mrlmer_full)
plot(mrlmer_full)
Anova(mrlmer_full)
qqnorm(resid(mrlmer_full))
qqline(resid(mrlmer_full)) 

# model without 2 way interaction
mrlmer_1 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc +
                   (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
Anova(mrlmer_1)

# model without myc_assoc
mrlmer_2 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc +
                   (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
summary(mrlmer_2)
Anova(mrlmer_2)

# model without incubation period
mrlmer_3 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate +
                   (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
Anova(mrlmer_3)
# model without isolate
mrlmer_4 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + 
                   (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
Anova(mrlmer_4)

# model without soil _moisture
mrlmer_5 <- lmer(lmt.mo ~ veg_type + pH + 
                   (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
Anova(mrlmer_5)

step(mrlmer_full)
# model chosen by step function
mrlmer_6 <- lmer(lmt.mo ~ veg_type + pH + isolate + incub_period.days +
                   (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
Anova(mrlmer_6)

# model comparison
anova(mrlmer_1,mrlmer_2,mrlmer_3,mrlmer_4,mrlmer_5,mrlmer_6)

# final model
mrlmer_fin <- lmer(lmt.mo ~ veg_type + pH + isolate + incub_period.days + 
                     (1 | incub_period.days:(isolate:exp)), data = fundecomp, REML=FALSE)
# visualizing the results 
summary(mrlmer_fin)
Anova(mrlmer_fin, type=2)
visreg(mrlmer_fin, "incub_period.days", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)
visreg(mrlmer_fin, "incub_period.days", by="veg_type",
       ylab='mass lost (proportion)', overlay=TRUE)
visreg(mrlmer_fin, "isolate", by="veg_type",
       ylab='mass lost (proportion)', overlay=TRUE)
visreg(mrlmer_fin, "pH", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)

# Decomposition differed the most in the prairie site. Isolates differed in their mass remaining with the exception of mel_white and Mort
# model diagnostics 
str(terms(mrlmer_fin))
isLMM(mrlmer_fin)
plot(mrlmer_fin)


# Just looking at the two isolates decomopsed at all of the sites (Mel_black and Mort)
mrlmer.mbm_full <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + isolate:incub_period.days + myc_assoc:incub_period.days + myc_assoc + soil_moisture:pH + myc_assoc:pH +
                      (1|incub_period.days:exp) ,data = fundecomp_melb.mort,  REML=FALSE)
summary(mrlmer.mbm_full)
Anova(mrlmer.mbm_full)
step(mrlmer.mbm_full)

mrlmer.mbm_step <- lmer(lmt.mo ~ veg_type + pH + isolate + incub_period.days +
                          (1|incub_period.days:exp) ,data = fundecomp_melb.mort,  REML=FALSE)
Anova(mrlmer.mbm_step)
mrlmer.mbm_mel <- lmer(lmt.mo ~ veg_type + pH + mel_status + incub_period.days +
                         (1|incub_period.days:exp),data = fundecomp_melb.mort,  REML=FALSE)
anova(mrlmer.mbm_full,mrlmer.mbm_step,mrlmer.mbm_mel)

# visualizing the results 

visreg(mrlmer.mbm_step, "incub_period.days", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)
visreg(mrlmer.mbm_step, "incub_period.days", by="veg_type",
       ylab='mass lost (proportion)', overlay=TRUE)
visreg(mrlmer.mbm_step, "isolate", by="veg_type",
       ylab='mass lost (proportion)', overlay=TRUE)
visreg(mrlmer.mbm_step, "pH", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)

