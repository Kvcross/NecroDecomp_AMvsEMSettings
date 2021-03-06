# title: Summer 2017 Fungal Necromass Collaboration Decay Constants
# author: KVB
#date: "9/26/2018"

# removing objects from the R Environment
rm(list=ls())

# load required packages
packload <- c("tidyverse","dplyr","broom", "ggplot2", "Rmisc", "broom", "MASS", "nlstools", "minpack.lm", "Hmisc", "car", "gtools", "lme4", "graphics", "rpart")
lapply(packload, library, character.only=TRUE)


# read in data file- with intitial values for MR
ifundecomp <- read.delim("data/NecroDecompinitial_summer2017.csv", sep = ",")

# removing NAs from mass remaining
ifundecomp <- ifundecomp %>% 
  filter(!is.na(mt.mo))

# subsetting
# Litter sp. decomposed at Cedar Creek in Praire fields
CCpr_ME <- filter(ifundecomp,  site  == "Prairie", isolate == "Mort")
CCpr_MB <- filter(ifundecomp,  site  == "Prairie", isolate == "Mel_Black")
CCpr_Ceno <- filter(ifundecomp,  site  == "Prairie", isolate == "Ceno")

# Litter sp. decomposed at Cedar Creek in Oak Savanna plots
CCsav_ME <- filter(ifundecomp,  site  == "Oak Savanna", isolate == "Mort")
CCsav_MB <- filter(ifundecomp,  site  == "Oak Savanna", isolate == "Mel_Black")
CCsav_MW <- filter(ifundecomp,  site  == "Oak Savanna", isolate == "Mel_White")  

# Litter sp. decomposed at Moores Creek in forest plots
MCfor_ME <- filter(ifundecomp,  site  == "Temperate Forest", isolate == "Mort")
MCfor_MB <- filter(ifundecomp,  site  == "Temperate Forest", isolate == "Mel_Black")

# calculating decay constants for each sp. in each ecosystem setting
# single Exponential Fit 
# k = ln(mt.mo)/(-time_year)
# mt.mo = the final mass/ the initial mass
# time_year = days incubated/ 365 days (units=yr-1)

# single Exponential Models for litter sp in Cedar Creek Prairie fields
# CC Pr ME
k.CCpr_ME = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=CCpr_ME)
summary(k.CCpr_ME) # parameter estimates and overall model fit
# CC Pr MB
k.CCpr_MB = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=CCpr_MB)
summary(k.CCpr_MB) # parameter estimates and overall model fit
# CC Pr Ceno
k.CCpr_Ceno = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=CCpr_Ceno)
summary(k.CCpr_Ceno) # parameter estimates and overall model fit

# plotting predicted values to look at line fit  
par(mfrow=c(1,3))
opar <- par(las = 1)
# CC Pr ME
tt <- seq(0, 6, length = 101)
plot(CCpr_ME$incub_period.wks, 
     CCpr_ME$mt.mo, main="Single Exp Model Fit CC Praire ME", 
     xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.CCpr_ME, list(incub_period.wks = tt)),col="red")
# CC Pr MB
tt <- seq(0, 6, length = 101)
plot(CCpr_MB$mt.mo ~ CCpr_MB$incub_period.wks,main="Single Exp Model Fit CC Praire MB",xlab="Time (weeks)", ylab="Mass Remaining",)
lines(tt,predict(k.CCpr_MB, list(incub_period.wks = tt)),col="red") +
text(cex=2)
# CC Pr Ceno
tt <- seq(0, 9, length = 101)
plot(CCpr_Ceno$mt.mo ~ CCpr_Ceno$incub_period.wks,main="Single Exp Model Fit CC Praire Ceno",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.CCpr_Ceno, list(incub_period.wks = tt)),col="red") +
text(cex=2)
par(opar)
dev.off()

# single Exponential Models for litter sp in Cedar Creek Oak Savanna Plots
# CC Oak Sav ME
k.CCsav_ME = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=CCsav_ME)
summary(k.CCsav_ME) # parameter estimates and overall model fit
# CC Oak Sav MB
k.CCsav_MB = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=CCsav_MB)
summary(k.CCsav_MB) # parameter estimates and overall model fit
# CC Oak Sav Ceno
k.CCsav_MW = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=CCsav_MW)
summary(k.CCsav_MW) # parameter estimates and overall model fit

# plotting predicted values to look at line fit  
par(mfrow=c(1,3))
opar <- par(las = 1)
# CC Oak Sav ME 
tt <- seq(0, 8, length = 101)
plot(CCsav_ME$incub_period.wks, CCsav_ME$mt.mo, main="Single Exp Model Fit CC Oak Sav ME",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.CCsav_ME, list(incub_period.wks = tt)),col="red") +
text(cex=2)
# CC Oak Sav MB black
tt <- seq(0, 8, length = 101)
plot(CCsav_MB$mt.mo ~ CCsav_MB$incub_period.wks,main="Single Exp Model Fit CC Oak Sav MB",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.CCsav_MB, list(incub_period.wks = tt)),col="red") +
text(cex=2)
# CC Oak Sav MB white
tt <- seq(0, 8, length = 101)
plot(CCsav_MW$mt.mo ~ CCsav_MW$incub_period.wks,main="Single Exp Model Fit CC Oak Sav MW",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.CCsav_MW, list(incub_period.wks = tt)),col="red")+
text(cex=2)
par(opar)
dev.off()

# single Exponential Models for litter sp in Cedar Creek Oak Savanna Plots
# MC Forest ME 
k.MCfor_ME = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=MCfor_ME)
summary(k.MCfor_ME) # parameter estimates and overall model fit
# CC MC Forest MB
k.MCfor_MB = nls(mt.mo~exp(-k*incub_period.wks), start=list(k=0.01), data=MCfor_MB)
summary(k.MCfor_MB) # parameter estimates and overall model fit

# plotting predicted values to look at line fit  
par(mfrow=c(1,2))
opar <- par(las = 1)
# MC forest ME 
tt <- seq(0, 12, length = 101)
plot(MCfor_ME$incub_period.wks, MCfor_ME$mt.mo, main="Single Exp Model Fit MC Forest ME",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.MCfor_ME, list(incub_period.wks = tt)),col="red")+
text(cex=2)
# CC forest MB
tt <- seq(0, 12, length = 101)
plot(MCfor_MB$mt.mo ~ MCfor_MB$incub_period.wks,main="Single Exp Model Fit MC Forest MB",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(k.MCfor_MB, list(incub_period.wks = tt)),col="red")+
text(cex=2)
# Double Exponential Fit
# mt.mo = S * e^-k1t + (1-S) * e^-k2t
# mt.mo = the final mass/ the initial mass
# k1 = k1 is the rate constant for the labile component
# k2 = k2 is the rate constant for the resistant component 
# t = days incubated/ 365 days (units=yr-1)
# s is the initial proportion of labile material
# 1-s is the initial proportion of resistant material 


# fitting a double exponential model for litter sp in Cedar Creek Praire fields
# using Levenberg-Marquardt algorithm to estimate starting values for nls model # using the minpak.lm package
# CC pr ME
guess.CCpr_ME = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=1,k1=1,k2=0.1), data=CCpr_ME) 
coef(guess.CCpr_ME)
K2.CCpr_ME = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.77,k1=1.06,k2=-0.03), data=CCpr_ME)
summary(K2.CCpr_ME) # parameter estimates and overall model fit
#CC pr MB
guess.CCpr_MB = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=-0.8,k1=0.3,k2=0.001), data=CCpr_MB) # check on this- not sure if I can fit a model to 3 time points
coef(guess.CCpr_MB)
K2.CCpr_MB = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=-1.7,k1=0.3,k2=0.001), data=CCpr_MB)
summary(K2.CCpr_MB) # parameter estimates and overall model fit
#CC pr Ceno
guess.CCpr_Ceno = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.5,k1=1,k2=1), data=CCpr_Ceno) 
coef(guess.CCpr_Ceno)
K2.CCpr_Ceno = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.76,k1=1.54,k2=0.041), data=CCpr_Ceno)
summary(K2.CCpr_Ceno) # parameter estimates and overall model fit

# plotting the model fits (double exp)
par(mfrow=c(1,3))
opar <- par(las = 1)
tt <- seq(0, 12, length = 101)
plot(CCpr_ME$incub_period.wks,CCpr_ME$mt.mo,main="Double Exp Model Fit CC Prairie ME",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.CCpr_ME, list(incub_period.wks = tt)),col="red")
plot(CCpr_MB$incub_period.wks,CCpr_MB$mt.mo,main="Double Exp Model Fit CC Prairie MB",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.CCpr_MB, list(incub_period.wks = tt)),col="red")
plot(CCpr_Ceno$incub_period.wks,CCpr_Ceno$mt.mo,main="Double Exp Model Fit CC Prairie  Ceno",xlab="Time (weeks)", ylab="Mass Remaining")
lines(CCpr_Ceno$incub_period.wks,predict(K2.CCpr_Ceno),col="red")
lines(tt,predict(K2.CCpr_Ceno, list(incub_period.wks = tt)),col="red")
par(opar)
dev.off()

# CC Oak Sav ME
guess.CCsav_ME = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.7,k1=1,k2=0.05), data=CCsav_ME) 
coef(guess.CCsav_ME)  
K2.CCsav_ME = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.2,k1=3.08,k2=0.09), data=CCsav_ME)
summary(K2.CCsav_ME)  # parameter estimates and overall model fit
#CC Oak Sav MB
guess.CCsav_MB = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.7,k1=1,k2=1), data=CCsav_MB) 
coef(guess.CCsav_MB)
K2.CCsav_MB = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.7,k1=0.7,k2=0.04), data=CCsav_MB)
summary(K2.CCsav_MB)  # parameter estimates and overall model fit
#CC Oak Sav MW
guess.CCsav_MW = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.7,k1=1,k2=1), data=CCsav_MW)
coef(guess.CCsav_MW)
K2.CCsav_MW = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.66,k1=1.6,k2=0.1), data=CCsav_MW)
summary(K2.CCsav_MW)  # parameter estimates and overall model fit

# plotting the model fits (double exp)
par(mfrow=c(1,3))
opar <- par(las = 1)
tt <- seq(0, 12, length = 101)
plot(CCsav_ME$incub_period.wks,CCsav_ME$mt.mo,main="Double Exp Model Fit CC Oak Sav ME",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.CCsav_ME, list(incub_period.wks = tt)),col="red")
plot(CCsav_MB$incub_period.wks,CCsav_MB$mt.mo,main="Double Exp Model Fit CC Oak Sav MB",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.CCsav_MB, list(incub_period.wks = tt)),col="red")
plot(CCsav_MW$incub_period.wks,CCsav_MW$mt.mo,main="Double Exp Model Fit CC Oak Sav MW",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.CCsav_MW, list(incub_period.wks = tt)),col="red")
par(opar)
dev.off()

#MC Forest ME
guess.MCfor_ME = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)),start=list(s=1,k1=1,k2=0.1),data=MCfor_ME) 
coef(guess.MCfor_ME)
K2.MCfor_ME = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.8,k1=1.2,k2=0.0001),data=MCfor_ME)
summary(K2.MCfor_ME)  # parameter estimates and overall model fit
#MC Forest MB
guess.MCfor_MB = nlsLM(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.7,k1=1,k2=0.1), data=MCfor_MB) 
coef(guess.MCfor_MB)
K2.MCfor_MB = nls(mt.mo~s*exp(-k1*incub_period.wks)+((1-s)*exp(-k2*incub_period.wks)), start=list(s=0.63,k1=0.627,k2=0.05), data=MCfor_MB)
summary(K2.MCfor_MB)  # parameter estimates and overall model fit

# plotting the model fits (double exp)
par(mfrow=c(1,2))
opar <- par(las = 1)
tt <- seq(0, 14, length = 101)
plot(MCfor_ME$incub_period.wks,MCfor_ME$mt.mo, main="Double Exp Model Fit MC Forest ME",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.MCfor_ME, list(incub_period.wks = tt)),col="red")
plot(MCfor_MB$incub_period.wks,MCfor_MB$mt.mo,main="Double Exp Model Fit MC Forest MB",xlab="Time (weeks)", ylab="Mass Remaining")
lines(tt,predict(K2.MCfor_MB, list(incub_period.wks = tt)),col="red")
par(opar)
# using AIC to compare model fits
aic.vals <- AIC(K2.CCsav_ME,K2.MCfor_ME,K2.CCsav_MB,K2.MCfor_MB)
write.table(aic.vals, "output/AIC_values.decaym.txt", sep = "\t", row.names = TRUE, quote = FALSE)

# using anova to compare single exponential model fits
anova(k.CCpr_ME,k.CCsav_ME,k.MCfor_ME,k.CCpr_MB,k.CCsav_MB,k.MCfor_MB)
anov.tblsing <- tidy(anov.singexp)
write.table(anov.tblsing, "output/anovatable.decaym.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# using anova to compare double exponential model fits
anov.doubexp <- anova(K2.CCsav_ME,K2.MCfor_ME,K2.CCsav_MB,K2.MCfor_MB)
anov.tbldoub <- tidy(anov.doubexp)
write.table(anov.tbldoub, "output/anovatable.decaym.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# looks like comparatively the double exponential models are better fits
# however, the std error looks pretty high for the estimates


# Putting model estimates into temporary data.frames:
model1Frame <- data.frame(Variable = rownames(summary(k.CCpr_ME)$coef),
                          Coefficient = summary(k.CCpr_ME)$coef[, 1],
                          SE = summary(k.CCpr_ME)$coef[, 2],
                          modelName = "Prairie Mortierella")
model2Frame <- data.frame(Variable = rownames(summary(K2.CCpr_ME)$coef),
                          Coefficient = summary(K2.CCpr_ME)$coef[, 1],
                          SE = summary(K2.CCpr_ME)$coef[, 2],
                          modelName = "Prairie Mortierella")
model3Frame <- data.frame(Variable = rownames(summary(k.CCsav_ME)$coef),
                          Coefficient = summary(k.CCsav_ME)$coef[, 1],
                          SE = summary(k.CCsav_ME)$coef[, 2],
                          modelName = "Oak Savanna Mortierella")
model4Frame <- data.frame(Variable = rownames(summary(K2.CCsav_ME)$coef),
                          Coefficient = summary(K2.CCsav_ME)$coef[, 1],
                          SE = summary(K2.CCsav_ME)$coef[, 2],
                          modelName = "Oak Savanna Mortierella")
model5Frame <- data.frame(Variable = rownames(summary(k.MCfor_ME)$coef),
                          Coefficient = summary(k.MCfor_ME)$coef[, 1],
                          SE = summary(k.MCfor_ME)$coef[, 2],
                          modelName = "Temperate Forest Mortierella")
model6Frame <- data.frame(Variable = rownames(summary(K2.MCfor_ME)$coef),
                          Coefficient = summary(K2.MCfor_ME)$coef[, 1],
                          SE = summary(K2.MCfor_ME)$coef[, 2],
                          modelName = "Temperate Forest Mortierella")
model7Frame <- data.frame(Variable = rownames(summary(k.CCpr_MB)$coef),
                          Coefficient = summary(k.CCpr_MB)$coef[, 1],
                          SE = summary(k.CCpr_MB)$coef[, 2],
                          modelName = "Prairie Meliniomyces Black")
model8Frame <- data.frame(Variable = rownames(summary(k.CCsav_MB)$coef),
                          Coefficient = summary(k.CCsav_MB)$coef[, 1],
                          SE = summary(k.CCsav_MB)$coef[, 2],
                          modelName = "Oak Savanna Meliniomyces Black")
model9Frame <- data.frame(Variable = rownames(summary(K2.CCsav_MB)$coef),
                          Coefficient = summary(K2.CCsav_MB)$coef[, 1],
                          SE = summary(K2.CCsav_MB)$coef[, 2],
                          modelName = "Oak Savanna Meliniomyces Black")
model10Frame <- data.frame(Variable = rownames(summary(k.MCfor_MB)$coef),
                          Coefficient = summary(k.MCfor_MB)$coef[, 1],
                          SE = summary(k.MCfor_MB)$coef[, 2],
                          modelName = "Temperate Forest Meliniomyces Black")
model11Frame <- data.frame(Variable = rownames(summary(K2.MCfor_MB)$coef),
                          Coefficient = summary(K2.MCfor_MB)$coef[, 1],
                          SE = summary(K2.MCfor_MB)$coef[, 2],
                          modelName = "Temperate Forest Meliniomyces Black")
model12Frame <- data.frame(Variable = rownames(summary(k.CCpr_Ceno)$coef),
                          Coefficient = summary(k.CCpr_Ceno)$coef[, 1],
                          SE = summary(k.CCpr_Ceno)$coef[, 2],
                          modelName = "Prairie Cenococcum")
model13Frame <- data.frame(Variable = rownames(summary(K2.CCpr_Ceno)$coef),
                           Coefficient = summary(K2.CCpr_Ceno)$coef[, 1],
                           SE = summary(K2.CCpr_Ceno)$coef[, 2],
                           modelName = "Prairie Cenococcum")
model14Frame <- data.frame(Variable = rownames(summary(k.CCsav_MW)$coef),
                           Coefficient = summary(k.CCsav_MW)$coef[, 1],
                           SE = summary(k.CCsav_MW)$coef[, 2],
                           modelName = "Oak Savanna Meliniomyces White")
model15Frame <- data.frame(Variable = rownames(summary(K2.CCsav_MW)$coef),
                           Coefficient = summary(K2.CCsav_MW)$coef[, 1],
                           SE = summary(K2.CCsav_MW)$coef[, 2],
                           modelName = "Oak Savanna Meliniomyces White")

# Combine these data frames into 1 data frame
combmodelframe <- data.frame(rbind(model1Frame, model2Frame, model3Frame,
                                   model4Frame, model5Frame, model6Frame,
                                   model7Frame, model8Frame, model9Frame,
                                   model10Frame, model11Frame, model12Frame,
                                   model13Frame, model14Frame, model15Frame)) 

combmodelframe.sing <- combmodelframe  %>%
  filter(!Variable %in% c("s")) 

# specifying the width of the95% confidence bands (The nls confidence intervals are 2.5% and 97.5% by default)
interval1 <- -qnorm((1-0.025)/2)  
interval2 <- -qnorm((1-0.97)/2)  

# plot
modcomp.plot <- ggplot(combmodelframe.sing, aes(colour = modelName)) +
       geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
       geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval1,
                                ymax = Coefficient + SE*interval1),
                            lwd = 1, position = position_dodge(width = 1/2)) +
       geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE") +
      coord_flip() +
      scale_x_discrete("Decay Constants") + 
      scale_y_continuous("Decay Constant Value") + 
      labs(color="Veg Community & Fungal Isolate") +
      ggtitle("Comparing Model Decay Constants (k)") +
      theme_classic(base_size=20) +
      theme(plot.title = element_text(hjust = 0.5))
ggsave(modcomp.plot, filename = "./figures/Decay_modelcomparison.png", units = "in", width = 15, height = 6, dpi = 300)


