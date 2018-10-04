# title: Summer 2017 Fungal Necromass Collaboration Mass Remaining Stats
# author: KVB
 #date: "9/21/2018"

# housekeeping 
rm(list=ls())

# load required packages
install.packages("pacman")
pacman::p_load(tidyverse,dplyr,broom, minpack.lm,Hmisc,car,gtools,lme4,nlme,minpack.lm)

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

fundecomp_melb.mort <- fundecomp  %>%
  filter(!isolate %in% c("Ceno","Mel_White")) 
  
# averages 
mass.avg <- fundecomp %>%
            select(veg_type,isolate, mass_remaining) %>%
            group_by(veg_type,isolate) %>%
            summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(mass.avg , "output/avgs_massdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# data exploration
str(fundecomp)
sum(is.na(fundecomp$mt.mo))

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

# doing a logistic transformation and then modelling linearly 
ggplot(fundecomp, aes(x = incub_period.wks, y = lmt.mo)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Logit(Mass Remaining") +
  xlab("Weeks of field Incubation)") +
  theme_classic() 

# looking for collinearity
vars <- c("soil_moisture", "pH", "incub_period.days","mt.mo")
pairs(fundecomp[,vars])
plot(fundecomp$soil_moisture ~ fundecomp$pH) 
cor.test(fundecomp$soil_moisture, fundecomp$pH) #low 0.3 collinearity between soil_moisture and pH

# using a linear mixed effect model to look at the effects of incubation period, myc association, veg_type and fungal isolate (fixed) and the effects of plot (random)  
# I will fit the most complex model and then eliminate fixed effects to determine if the eliminated effect is significant
# full model
mrlmer_full <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc + veg_type:isolate + veg_type:incub_period.days + incub_period.days:isolate + incub_period.days:isolate:veg_type +(1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_full)
Anova(mrlmer_full)
plot(mrlmer_full)
qqnorm(resid(mrlmer_full))
qqline(resid(mrlmer_full)) 

# model without 3 way interaction
mrlmer_r1 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc + veg_type:isolate + veg_type:incub_period.days + incub_period.days:isolate + (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r1)

# model without 2 way interaction incubation period:isolate
mrlmer_r2 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc + veg_type:isolate + veg_type:incub_period.days +  (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r2)
Anova(mrlmer_r2)

# model without 2 way interaction veg_type:incubation period
mrlmer_r3 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc + veg_type:isolate + (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r3)
Anova(mrlmer_r3)

# model without 2 way interaction veg_type:isolate
mrlmer_r4 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc + (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r4)
Anova(mrlmer_r4)

# model without myc association 
mrlmer_r5 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r5)
Anova(mrlmer_r5)

# model without incub_period.days
mrlmer_r6 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate +  (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r6)
Anova(mrlmer_r6)

# model without isolate
mrlmer_r7 <- lmer(lmt.mo ~ veg_type + pH + soil_moisture + (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r7)
Anova(mrlmer_r7)

# model without soil moisture
mrlmer_r8 <- lmer(lmt.mo ~ veg_type + pH +(1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r8)
Anova(mrlmer_r8)

# model without pH
mrlmer_r9 <-lmer(lmt.mo ~ veg_type + (1|unique_plotID), data = fundecomp,REML=FALSE)
summary(mrlmer_r9)
Anova(mrlmer_r9)

# model comparison
anova(mrlmer_full,mrlmer_r1,mrlmer_r2,mrlmer_r3,mrlmer_r4,mrlmer_r5,mrlmer_r6,mrlmer_r7,mrlmer_r8)

# final model
mrlmer_fin <-  lmer(lmt.mo ~ veg_type + pH + soil_moisture + isolate + incub_period.days + myc_assoc + veg_type:isolate + veg_type:incub_period.days + incub_period.days:isolate + (1|unique_plotID), data = fundecomp)
summary(mrlmer_fin)
mrlmer_fin <- Anova(mrlmer_fin, type=3)
mixmod <- tidy(mrlmer_fin)
sig <- tidy(mrlmer_fin)
write.table(mixmod, "output/mixmod.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sig, "output/mixmod.sig.txt", sep = "\t", row.names = FALSE, quote = FALSE)

