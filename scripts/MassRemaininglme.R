# title: Summer 2017 Fungal Necromass Collaboration Mass Remaining Analysis
# author: KVB
 #date: "9/21/2018"

# housekeeping 
rm(list=ls())

# load required packages
packload <- c("tidyverse","dplyr","broom","lmerTest","minpack.lm","Hmisc","car","gtools","lme4","nlme","minpack.lm","visreg", "lsmeans", "survival")
lapply(packload, library, character.only=TRUE)

# read in data file
fundecomp <- read.delim("data/CompiledNecroDecomp_summer2017.csv", sep = ",")

# specifying factors for data file
fundecomp  <-  within(fundecomp , {
            location = factor(location, ordered=FALSE)
            site = factor(site, ordered=FALSE)
            exp = factor(exp, ordered=FALSE)
            location = factor(location, ordered=FALSE)
            plot_Id_pair = factor(plot_Id_pair, ordered=FALSE)
            unique_plotID = factor(unique_plotID, ordered=FALSE)
            treatment = factor(treatment, ordered=FALSE)
            myc_assoc  = factor(myc_assoc, ordered=FALSE)
            isolate = factor(isolate, ordered=FALSE)
            incub_period.days = factor(incub_period.days, ordered=FALSE)
            incub_period.wks = factor(incub_period.wks , ordered=FALSE)
            incub_period.months = factor(incub_period.months , ordered=FALSE)
            incub_period.years = factor(incub_period.years, ordered=FALSE)
})

# removing NAs for mass remaining
fundecomp <- fundecomp %>% 
            filter(!is.na(mt.mo))
sum(is.na(fundecomp$mt.mo))
  
# data exploration
str(fundecomp)
summary(fundecomp)
# count the data to confirm that the design is unbalanced 
ftable(xtabs(~ isolate + incub_period.days + site, fundecomp))

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
  filter(location != "moores creek") # just the cedar creek experiments 

MCfd_melb.mort <- fd_melb.mort %>%
  filter(location != "cedar creek") # Just the moores creek experiment/ temperate forest veg type

CCfd_os <- CCfd_melb.mort %>%
  filter(site != "Prairie")  # Just the oak savanna veg type

CCfd_p <- CCfd_melb.mort %>%
  filter(site != "Oak Savanna") # Just the prairie veg type 

CCfd_pex2 <- CCfd_p %>%
  filter(!exp %in% "CC_praire1")
write.csv(CCfd_pex2, file = "CCPrex2.csv")

CCfd_pf5 <- CCfd_melb.mort %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(!myc_assoc %in% "mixed") %>% 
  filter(!location.1%in% c("field_0","field_80"))

MREarlyd <- fd_melb.mort  %>%
  filter(!exp %in% "CC_praire1") %>%
  filter(incub_period.days %in% c("11","14"))
         
# averages for all isolates
mass.avg1 <- fd_melb.mort %>%
    select(site,isolate,incub_period.days, mass_remaining) %>%
    group_by(site,isolate,incub_period.days) %>%
    summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(mass.avg1 , "output/avgs_massdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

mass.avg2 <- fd_melb.mort %>%
  select(site, myc_assoc, mass_remaining) %>%
  group_by(site, myc_assoc) %>%
  summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(mass.avg2 , "output/avgs_massdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# doing a logistic transformation and then modelling linearly using 3 seperate linear mixed effects model using plot (CC) or plot pairs (MC) as random effects
# for the remainder of the analyis I will exclude Ceno and Mel White/ The mixed treatment from the Oak Savanna vegetation type (datasheet = fundecomp_melb.mort_myc)
# lsmeans will be used for pairwise comparisons 
# to compare isolates across sites I will generate K values-refer to DecayConstants.R script

# data visualization
ggplot(fd_melb.mort, aes(x = incub_period.days, y = lmt.mo, color=isolate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Logit(Mass Remaining)") +
  xlab("Days of field Incubation)") +
  facet_grid(. ~ site) +
  theme_classic()

# looking for collinearity
vars <- c("soil_moisture", "pH", "incub_period.days","mt.mo")
pairs(fd_melb.mort[,vars])
plot(fd_melb.mort$soil_moisture ~ fd_melb.mortc$pH) 
cor.test(fd_melb.mort$soil_moisture, fd_melb.mort$pH) #low 0.42 collinearity between soil_moisture and pH

### Moores Creek lme

mclmer <- lmer(lmt.mo ~  myc_assoc * isolate * incub_period.wks + soil_moisture + (1|plot_Id_pair),
               data = MCfd_melb.mort,na.action = na.exclude)
summary(mclmer)
# looking at residuals
plot(mclmer)
qqnorm(resid(mclmer))
qqline(resid(mclmer)) 
hist(resid(mclmer))
# anova results for fixed effects
MC_f_lme <- anova(mclmer,type=3,ddf="Kenward-Roger") 
# using the Kenward-Roger denominator degrees of freedom approximation which is more conservitive for unbalanced mixed model 
# anova like results for random effect
MC_f_rlme <- ranova(mclmer) 

# visualizing the results
visreg(mclmer, "incub_period.wks", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)

visreg(mclmer, "isolate", by="myc_assoc",
       ylab='mass lost (proportion)', overlay=TRUE)

# Moores Creek ls means comparisons myc assoc
mc.myc = lsmeans(mclmer, ~ myc_assoc)
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

# Moores Creek ls means comparisons isolate
mc.iso = lsmeans(mclmer, ~ isolate)
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

# ls means comparisons incubation time 
mc.tim = lsmeans(mclmer, ~ incub_period.wks)
pairs(mc.tim,adjust="tukey")
mtimpw = cld(mc.tim,alpha = 0.05,Letters = letters, adjust  = "tukey") 
mtimpw 

# Plot
ggplot(mtimpw ,aes(x = incub_period.wks, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# ls means comparisons isolate : incubation time
mc.timis = lsmeans(mclmer, ~ isolate + incub_period.wks)
pairs(mc.timis,adjust="tukey")
mtimispw = cld(mc.timis,alpha = 0.05,Letters = letters, adjust  = "tukey") 
mtimispw 

# Plot
ggplot(mtimispw ,aes(x = incub_period.wks, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# Prairie 
# Only using data from Experiment 2- could use data from same field across the two experiments 
sum(is.na(CCfd_p2$soil_moisture)) # 22 missing entries out of 47 for soil moisture so I am leaving it out of the model

# Cedar Creek - prairie lme
ccplmer <- lmer(lmt.mo ~ myc_assoc * isolate * incub_period.wks + (1|plot_Id_pair),
               data = CCfd_pex2, na.action = na.exclude)
summary(ccplmer)
plot(ccplmer)
qqnorm(resid(ccplmer))
qqline(resid(ccplmer)) 
hist(resid(ccplmer))
# anova results for fixed effects
CC_p_lme <- anova(ccplmer,type=3, ddf="Kenward-Roger")
# anova like results for random effect
CC_p_rlme <- ranova(ccplmer) 

# visualizing the results
visreg(ccplmer, "incub_period.wks", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)

visreg(ccplmer, "isolate", by="myc_assoc",
       ylab='mass lost (proportion)', overlay=TRUE)

# ls means comparisons myc assoc
cp.myc = lsmeans(ccplmer, ~ myc_assoc)
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

# ls means comparisons isolate
cp.iso = lsmeans(ccplmer, ~ isolate)
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
cp.tim = lsmeans(ccplmer, ~ incub_period.wks +isolate)
pairs(cp.tim,adjust="tukey")
cptimpw = cld(cp.tim,alpha = 0.05,Letters = letters, adjust  = "tukey") 
cptimpw 

# Plot
ggplot(cptimpw ,aes(x = incub_period.wks, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# Oak Savanna
# Cedar Creek oak savanna lme
ccoslmer <- lmer(lmt.mo ~ myc_assoc * isolate * incub_period.wks + soil_moisture  + (1|plot_Id_pair),
                data = CCfd_os, na.action = na.exclude)
summary(ccoslmer)
plot(ccoslmer)
qqnorm(resid(ccoslmer))
qqline(resid(ccoslmer)) 
hist(resid(ccoslmer))
# anova results for fixed effects
CC_os_lme <-  anova(ccoslmer,type=3, ddf="Kenward-Roger")
# anova like results for random effect
CC_os_rlme <- ranova(ccoslmer) 

# visualizing the results
visreg(ccoslmer, "incub_period.wks", by="isolate",
       ylab='mass lost (proportion)', overlay=TRUE)

visreg(ccoslmer, "isolate", by="myc_assoc",
       ylab='mass lost (proportion)', overlay=TRUE)

# ls means comparisons myc assoc
cp.myc = lsmeans(ccoslmer, ~ myc_assoc)
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

# ls means comparisons isolate
cos.iso = lsmeans(ccoslmer, ~ isolate)
pairs(cos.iso,adjust="tukey")
cosisopw = cld(cp.iso,alpha = 0.05,Letters = letters, adjust  = "tukey") 
cosisopw 

# Comparing early stages of decay
MREd.mod = lm(lmt.mo ~  site * myc_assoc * isolate + soil_moisture, data=MREarlyd, na.action=na.exclude)
summary(MREd.mod)
plot(MREd.mod) 
Anova(MREd.mod,type = "III")
qqnorm(resid(MREd.mod))
qqline(resid(MREd.mod))

# ls means comparisons site : isolate
early.comp = lsmeans(MREd.mod, ~ site + isolate)
pairs(early.comp,adjust="tukey")
early.comppw = cld(early.comp,alpha = 0.05,Letters = letters, adjust  = "tukey") 
early.comppw 

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
cos.tim = lsmeans(ccoslmer, ~ incub_period.wks)
pairs(cos.tim,adjust="tukey")
costimpw = cld(cos.tim,alpha = 0.05,Letters = letters, adjust  = "tukey") 
costimpw 

# Plot
ggplot(costimpw ,aes(x = incub_period.wks, y = lsmean, color = .group)) +
  geom_point(shape  = 15, size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")


# ls means comparisons isolate : incubation time
cos.timis = lsmeans(ccoslmer, ~ isolate + incub_period.wks)
pairs(cos.timis,adjust="tukey")
costimispw = cld(cos.timis,alpha = 0.05,Letters = letters, adjust  = "tukey") 
costimispw 

# Plot
ggplot(costimispw ,aes(x = incub_period.wks, y = lsmean, color = .group, shape=isolate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n mass remaining (logit transformed)")

# saving lmer models as one table
tidy_tf <- tidy(MC_f_lme)
tidy_pr <- tidy(CC_p_lme)
tidy_os <- tidy(CC_os_lme)
tidy_tf_r <- tidy(MC_f_rlme)
tidy_pr_r <- tidy(CC_p_rlme)
tidy_os_r <- tidy(CC_os_rlme)

write.table(tidy_tf, "output/forestLME.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tidy_pr, "output/prairieLME.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tidy_os, "output/oaksavLME.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tidy_tf_r, "output/forestRLME.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tidy_pr_r, "output/prairieRLME.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tidy_os_r, "output/oaksavRLME.txt", sep = "\t", row.names = FALSE, quote = FALSE)




