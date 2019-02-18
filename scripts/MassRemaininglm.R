# title: Summer 2017 Fungal Necromass Collaboration Mass Remaining Analysis
# author: KVB
 #date: "9/21/2018"

# housekeeping 
rm(list=ls())

# load required packages
packload <- c("tidyverse","dplyr","broom","lmerTest","minpack.lm","Hmisc","car","gtools","lme4","nlme","minpack.lm","visreg", "lsmeans")
lapply(packload, library, character.only=TRUE)

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

# removing NAs for mass remaining
fundecomp <- fundecomp %>% 
            filter(!is.na(mt.mo))

sum(is.na(fundecomp$mt.mo))

  
# averages for all isolates
mass.avg <- fundecomp %>%
            select(veg_type,isolate, mass_remaining) %>%
            group_by(veg_type,isolate) %>%
            summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(mass.avg , "output/avgs_massdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# data exploration
str(fundecomp)
summary(fundecomp)
# count the data to confirm that the design is unbalanced 
ftable(xtabs(~ isolate + incub_period.days + veg_type, fundecomp))

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

# subsetting
fd_myc <- fundecomp  %>%
  filter(myc_assoc != "mixed") # removing the mixed mycorrhizal association treatment 

fd_melb.mort <- fd_myc  %>%
  filter(isolate !="Ceno" & isolate !="Mel_White") # removing the Ceno and Mel_white isolates

CCfd_melb.mort <- fd_melb.mort %>%
  filter(site != "moores creek") # just the cedar creek experiments 

MCfd_melb.mort <- fd_melb.mort %>%
  filter(site != "cedar creek") # Just the moores creek experiment/ temperate forest veg type

CCfd_os <- CCfd_melb.mort %>%
  filter(veg_type != "Prairie")  # Just the oak savanna veg type

CCfd_p <- CCfd_melb.mort  %>%
  filter(veg_type != "Oak Savanna") # Just the prairie veg type 

# doing a logistic transformation and then modelling linearly (I tried using a linear mixed effects model, but I don't have a truly random factor for the 3 combined experiments)
# for the remainder of the analyis I will exclude Ceno and Mel White/ The mixed treatment from the Oak Savanna vegetation type (datasheet = fundecomp_melb.mort_myc)
# I will split the analysis up for the different experiments and sites (to account for differences in climate between study sites in IN & MN) 
# and then try a model where all of the data is included
# lsmeans will be used for pairwise comparisons 
# to compare isolates across sites I will generate K values-refer to DecayConstants.R script

# data visualization
ggplot(fd_melb.mort, aes(x = incub_period.days, y = lmt.mo, color=isolate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Logit(Mass Remaining)") +
  xlab("Days of field Incubation)") +
  facet_grid(. ~ veg_type) +
  theme_classic()

# looking for collinearity
vars <- c("soil_moisture", "pH", "incub_period.days","mt.mo")
pairs(fd_melb.mort[,vars])
plot(fd_melb.mort$soil_moisture ~ fd_melb.mort$pH) 
cor.test(fd_melb.mort$soil_moisture, fd_melb.mort$pH) #low 0.42 collinearity between soil_moisture and pH

# Linear models - Separating Moores Creek and Cedar Creek  Experiments 
# moores creek 
mcmrlm_full <- lm(lmt.mo ~ myc_assoc +  isolate + factor(incub_period.days) + pH + soil_moisture + 
                    isolate:factor(incub_period.days) +  isolate: myc_assoc + soil_moisture:pH + myc_assoc:pH + myc_assoc:soil_moisture,
                  data = MCfd_melb.mort,  na.action=na.exclude)
summary(mcmrlm_full)
plot(mcmrlm_full) 
Anova(mcmrlm_full,type = "II")
qqnorm(resid(mcmrlm_full))
qqline(resid(mcmrlm_full)) 

# I am missing pH for the 14 day incubation period
# removing pH /interactions that are not significant
mcmrlm <- lm(lmt.mo ~  myc_assoc + soil_moisture + isolate + factor(incub_period.days) + isolate:incub_period.days,
             data = MCfd_melb.mort,  na.action=na.exclude)
summary(mcmrlm)
plot(mcmrlm_full) 
Anova(mcmrlm,type = "II")
qqnorm(resid(mcmrlm))
qqline(resid(mcmrlm)) 

# Plot
# ls means comparisons myc assoc
mc.myc = lsmeans(mcmrlm, ~ myc_assoc)
pairs(mc.myc,adjust="tukey")
mmycpw = cld(mc.myc,alpha = 0.05,Letters = letters, adjust  = "tukey") 
mmycpw 

ggplot(mmycpw,aes(x = myc_assoc, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons isolate
mc.iso = lsmeans(mcmrlm, ~ isolate)
pairs(mc.iso,adjust="tukey")
misopw = cld(mc.iso,alpha = 0.05,Letters = letters, adjust  = "tukey") 
misopw 

# Plot
ggplot(misopw ,aes(x = isolate, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons incubation time & isolate
mc.tim = lsmeans(mcmrlm, ~ incub_period.days)
pairs(mc.tim,adjust="tukey")
mtimpw = cld(mc.tim,alpha = 0.05,Letters = letters, adjust  = "tukey") 
mtimpw 

# Plot
ggplot(mtimpw ,aes(x = incub_period.days, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# cedar creek-prairie

sum(is.na(CCfd_p$pH))
# pH and soil moisture measurments are missing for 46 of the 92 mass remaining entries 
# leaving pH and soil moisture out of the model
cpmrlm_full <- lm(lmt.mo ~ myc_assoc + location + isolate + factor(incub_period.days) + isolate:factor(incub_period.days) + isolate: myc_assoc + location:isolate,
                   data = CCfd_p,  na.action=na.exclude)
summary(cpmrlm_full)
plot(cpmrlm_full) 
Anova(cpmrlm_full,type = "II")
qqnorm(resid(cpmrlm_full))
qqline(resid(cpmrlm_full)) 

# removing interactions that are not significant
cpmrlm <- lm(lmt.mo ~ myc_assoc + location + isolate + factor(incub_period.days),
             data = CCfd_p,  na.action=na.exclude)
summary(cpmrlm)
plot(cpmrlm)
Anova(cpmrlm,type = "II")
qqnorm(resid(cpmrlm))
qqline(resid(cpmrlm)) 

# ls means comparisons myc assoc
cp.myc = lsmeans(cpmrlm, ~ myc_assoc)
pairs(cp.myc,adjust="tukey")
pmycpw = cld(cp.myc,alpha = 0.05,Letters = letters, adjust  = "tukey") 
pmycpw 

ggplot(pmycpw,aes(x = myc_assoc, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons location
cp.loc = lsmeans(cpmrlm, ~ location)
pairs(cp.loc,adjust="tukey")
cplocpw = cld(cp.loc,alpha = 0.05,Letters = letters, adjust  = "tukey") 
cplocpw 

# Plot
ggplot(cplocpw ,aes(x = location, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons isolate
cp.iso = lsmeans(cpmrlm, ~ isolate)
pairs(cp.iso,adjust="tukey")
cpisopw = cld(cp.iso,alpha = 0.05,Letters = letters, adjust  = "tukey") 
cpisopw 

# Plot
ggplot(cpisopw,aes(x = isolate, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons incubation time 
cp.tim = lsmeans(cpmrlm, ~ incub_period.days)
pairs(cp.tim,adjust="tukey")
cptimpw = cld(cp.tim,alpha = 0.05,Letters = letters, adjust  = "tukey") 
cptimpw 

# Plot
ggplot(cptimpw ,aes(x = incub_period.days, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# cedar creek-Oak Savanna

cosmrlm_full <- lm(lmt.mo ~ location +  isolate + factor(incub_period.days) + pH + soil_moisture + 
                    isolate:factor(incub_period.days) +  isolate:location + soil_moisture:pH + location:pH + location:soil_moisture,
                  data = CCfd_os,  na.action=na.exclude)
summary(cosmrlm_full)
plot(cosmrlm_full) 
Anova(cosmrlm_full,type = "II")
qqnorm(resid(cosmrlm_full))
qqline(resid(cosmrlm_full)) 

# Removing non-significant interactions

cosmrlm <- lm(lmt.mo ~ location +  isolate + factor(incub_period.days) + pH + soil_moisture + 
                isolate:factor(incub_period.days) +  isolate:location + soil_moisture:pH + location:pH + location:soil_moisture,
              data = CCfd_os,  na.action=na.exclude)

summary(cosmrlm_full)
plot(cosmrlm_full) 
Anova(cosmrlm_full,type = "II")
qqnorm(resid(cosmrlm_full))
qqline(resid(cosmrlm_full)) 

# removing interactions that are not significant
cosmrlm <- lm(lmt.mo ~ location +  isolate + factor(incub_period.days) + pH + soil_moisture,
              data = CCfd_os,  na.action=na.exclude)
summary(cosmrlm)
plot(cosmrlm)
Anova(cosmrlm,type = "II")
qqnorm(resid(cosmrlm))
qqline(resid(cosmrlm)) 


# ls means comparisons location
cos.loc = lsmeans(cosmrlm, ~ location)
pairs(cos.loc,adjust="tukey")
coslocpw = cld(cos.loc,alpha = 0.05,Letters = letters, adjust  = "tukey") 
coslocpw 

# Plot
ggplot(coslocpw ,aes(x = location, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons isolate
cos.iso = lsmeans(cosmrlm, ~ isolate)
pairs(cos.iso,adjust="tukey")
cosisopw = cld(cp.iso,alpha = 0.05,Letters = letters, adjust  = "tukey") 
cosisopw 

# Plot
ggplot(cosisopw,aes(x = isolate, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons incubation time 
cos.tim = lsmeans(cosmrlm, ~ incub_period.days)
pairs(cos.tim,adjust="tukey")
costimpw = cld(cos.tim,alpha = 0.05,Letters = letters, adjust  = "tukey") 
costimpw 

# Plot
ggplot(costimpw ,aes(x = incub_period.days, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")







