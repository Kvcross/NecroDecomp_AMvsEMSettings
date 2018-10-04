# title: Summer 2017 Fungal Necromass Collaboration Decay Constants
# author: KVB
#date: "9/26/2018"

# removing objects from the R Environment
rm(list=ls())

# load required packages
install.packages("pacman")
pacman::p_load(tidyverse,dplyr,broom, minpack.lm,Hmisc,car,gtools,lme4,nlme,minpack.lm)

# read in data file- with intitial values for MR
ifundecomp <- read.delim("data/NecroDecompinitial_summer2017.csv", sep = ",")

# removing NAs from mass remaining
ifundecomp <- ifundecomp %>% 
  filter(!is.na(mt.mo))

# subsetting
# Litter sp. decomposed at Cedar Creek in Praire fields
CCpr_ME <- filter(ifundecomp,  veg_type  == "Prairie", isolate == "Mort")
CCpr_MB <- filter(ifundecomp,  veg_type  == "Prairie", isolate == "Mel_Black")
CCpr_Ceno <- filter(ifundecomp,  veg_type  == "Prairie", isolate == "Ceno")

# Litter sp. decomposed at Moores Creek in forest plots
MCfor_ME <- filter(ifundecomp,  veg_type  == "Forest", isolate == "Mort")
MCfor_MB <- filter(ifundecomp,  veg_type  == "Forest", isolate == "Mel_Black")

# Litter sp. decomposed at Cedar Creek in Oak Savanna plots
CCsav_ME <- filter(ifundecomp,  veg_type  == "Oak Savanna", isolate == "Mort")
CCsav_MB <- filter(ifundecomp,  veg_type  == "Oak Savanna", isolate == "Mel_Black")
CCsav_MW <- filter(ifundecomp,  veg_type  == "Oak Savanna", isolate == "Mel_White")  

# calculating decay constants for each sp. in each ecosystem setting
# single Exponential Fit 
# k = ln(mt.mo)/(-time_year)
# mt.mo = the final mass/ the initial mass
# time_year = days incubated/ 365 days (units=yr-1)

# single Exponential Models for litter sp in Cedar Creek Prairie fields
# CC Pr ME
k.CCpr_ME = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=CCpr_ME)
summary(k.CCpr_ME)
# CC Pr MB
k.CCpr_MB = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=CCpr_MB)
summary(k.CCpr_MB)
# CC Pr Ceno
k.CCpr_Ceno = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=CCpr_Ceno)
summary(k.CCpr_Ceno)

# plotting predicted values to look at line fit  
par(mfrow=c(1,3))
# CC Pr ME
plot(CCpr_ME$incub_period.years, CCpr_ME$mt.mo, main="Single Exp Model Fit CC Praire ME")
lines(CCpr_ME$incub_period.years,predict(k.CCpr_ME),col="red")
# CC Pr MB
plot(CCpr_MB$mt.mo ~ CCpr_MB$incub_period.years,main="Single Exp Model Fit CC Praire MB")
lines(CCpr_MB$incub_period.years,predict(k.CCpr_MB),col="red")
# CC Pr Ceno
plot(CCpr_Ceno$mt.mo ~ CCpr_Ceno$incub_period.years,main="Single Exp Model Fit CC Praire Ceno")
lines(CCpr_Ceno$incub_period.years,predict(k.CCpr_Ceno),col="red")
dev.off()

# single Exponential Models for litter sp in Cedar Creek Oak Savanna Plots
# CC Oak Sav ME
k.CCsav_ME = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=CCsav_ME)
summary(k.CCsav_ME)
# CC Oak Sav MB
k.CCsav_MB = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=CCsav_MB)
summary(k.CCsav_MB)
# CC Oak Sav Ceno
k.CCsav_MW = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=CCsav_MW)
summary(k.CCsav_MW)

# plotting predicted values to look at line fit  
par(mfrow=c(1,3))
# CC Oak Sav ME 
plot(CCsav_ME$incub_period.years, CCsav_ME$mt.mo, main="Single Exp Model Fit CC Oak Sav ME")
lines(CCsav_ME$incub_period.years,predict(k.CCsav_ME),col="red")
# CC Oak Sav MB
plot(CCsav_MB$mt.mo ~ CCsav_MB$incub_period.years,main="Single Exp Model Fit CC Oak Sav MB")
lines(CCsav_MB$incub_period.years,predict(k.CCsav_MB),col="red")
# CC Oak Sav Ceno
plot(CCsav_MW$mt.mo ~ CCsav_MW$incub_period.years,main="Single Exp Model Fit CC Oak Sav MW")
lines(CCsav_MW$incub_period.years,predict(k.CCsav_MW),col="red")
dev.off()

# single Exponential Models for litter sp in Cedar Creek Oak Savanna Plots
# MC Forest ME 
k.MCfor_ME = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=MCfor_ME)
summary(k.MCfor_ME)
# CC MC Forest MB
k.MCfor_MB = nls(mt.mo~exp(-k*incub_period.years), start=list(k=0.01), data=MCfor_MB)
summary(k.MCfor_MB)

# plotting predicted values to look at line fit  
par(mfrow=c(1,2))
# MC forest ME 
plot(MCfor_ME$incub_period.years, MCfor_ME$mt.mo, main="Single Exp Model Fit MC Forest ME")
lines(MCfor_ME$incub_period.years,predict(k.MCfor_ME),col="red")
# CC forest MB
plot(MCfor_MB$mt.mo ~ MCfor_MB$incub_period.years,main="Single Exp Model Fit MC Forest MB")
lines(MCfor_MB$incub_period.years,predict(k.MCfor_MB),col="red")


# fitting a double exponential model for litter sp in Cedar Creek Praire fields
# using Levenberg-Marquardt algorithm to estimate starting values for nls model # using the minpak.lm package
# CC pr ME
guess.CCpr_ME = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=CCpr_ME) # check on this
coef(guess.CCpr_ME)
K2.CCpr_ME = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.77,k=55,h=-1.5), data=CCpr_ME)
summary(K2.CCpr_ME)
#CC pr MB
guess.CCpr_MB = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=-0.8,k=17,h=17), data=CCpr_MB) # check on this
coef(guess.CCpr_MB)
K2.CCpr_MB = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=-0.8,k=17,h=17), data=CCpr_MB)
summary(K2.CCpr_MB)
#CC pr Ceno
guess.CCpr_Ceno = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=CCpr_Ceno) # check on this
coef(guess.CCpr_Ceno)
K2.CCpr_Ceno = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.75,k=80.7,h=2.1), data=CCpr_Ceno)
summary(K2.CCpr_Ceno)

par(mfrow=c(1,3))
plot(CCpr_ME$incub_period.years,CCpr_ME$mt.mo,main="Double Exp Model Fit CC Prairie ME")
lines(CCpr_ME$incub_period.years,predict(K2.CCpr_ME),col="red")
plot(CCpr_MB$incub_period.years,CCpr_MB$mt.mo,main="Double Exp Model Fit CC Prairie MB")
lines(CCpr_MB$incub_period.years,predict(guess.CCpr_MB),col="red")
plot(CCpr_Ceno$incub_period.years,CCpr_Ceno$mt.mo,main="Double Exp Model Fit CC Prairie  MW")
lines(CCpr_Ceno$incub_period.years,predict(K2.CCpr_Ceno),col="red")
dev.off()

# CC Oak Sav ME
guess.CCsav_ME = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=CCsav_ME) # check on this- not sure if I can fit a model to 3 time points
coef(guess.CCsav_ME)
K2.CCsav_ME = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.70,k=160,h=4.6), data=CCsav_ME)
summary(K2.CCsav_ME)
#CC Oak Sav MB
guess.CCsav_MB = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=CCsav_MB) # check on this
coef(guess.CCsav_MB)
K2.CCsav_MB = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.71,k=36,h=2), data=CCsav_MB)
summary(K2.CCsav_MB)
#CC Oak Sav MW
guess.CCsav_MW = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=CCsav_MW) # check on this
coef(guess.CCsav_MW)
K2.CCsav_MW = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.65,k=84.7,h=5.5), data=CCsav_MW)
summary(K2.CCsav_MW)

par(mfrow=c(1,3))
plot(CCsav_ME$incub_period.years,CCsav_ME$mt.mo,main="Double Exp Model Fit CC Oak Sav ME")
lines(CCsav_ME$incub_period.years,predict(K2.CCsav_ME),col="red")
plot(CCsav_MB$incub_period.years,CCsav_MB$mt.mo,main="Double Exp Model Fit CC Oak Sav MB")
lines(CCsav_MB$incub_period.years,predict(K2.CCsav_MB),col="red")
plot(CCsav_MW$incub_period.years,CCsav_MW$mt.mo,main="Double Exp Model Fit CC Oak Sav MW")
lines(CCsav_MW$incub_period.years,predict(K2.CCsav_MW),col="red")
dev.off()


#MC Forest ME
guess.MCfor_ME = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=MCfor_ME) 
coef(guess.MCfor_ME)
K2.MCfor_ME = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.84,k=48,h=0.03), data=MCfor_ME)
summary(K2.MCfor_ME)
#MC Forest MB
guess.MCfor_MB = nlsLM(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=1,k=1,h=1), data=MCfor_MB) # check on this
coef(guess.MCfor_MB)
K2.MCfor_MB = nls(mt.mo~s*exp(-k*incub_period.years)+((1-s)*exp(-h*incub_period.years)), start=list(s=0.63,k=23.7,h=1.8), data=MCfor_MB)
summary(K2.MCfor_MB)

# plotting the output
par(mfrow=c(1,2))
plot(MCfor_ME$incub_period.years,MCfor_ME$mt.mo, main="Double Exp Model Fit MC Forest ME")
lines(MCfor_ME$incub_period.years,predict(K2.MCfor_ME),col="red")
plot(MCfor_MB$incub_period.years,MCfor_MB$mt.mo,main="Double Exp Model Fit MC Forest MB")
lines(MCfor_MB$incub_period.years,predict(K2.MCfor_MB),col="red")

# using AIC to compare model fits
AIC(k.CCpr_ME,K2.CCpr_ME,k.CCsav_ME,K2.CCsav_ME,k.MCfor_ME,K2.MCfor_ME,
    k.CCpr_MB,guess.MCfor_MB,k.CCsav_MB,K2.CCsav_MB,k.MCfor_MB,K2.MCfor_MB,
    k.CCpr_Ceno,K2.CCpr_Ceno,k.CCsav_MW,K2.CCsav_MW) 

# looks like comparatively the double exponential models are better fits- now to compare across veg types
AIC(K2.CCpr_ME,K2.CCsav_ME,K2.MCfor_ME,K2.CCsav_MW,K2.CCsav_MB,K2.MCfor_MB,K2.CCpr_Ceno) 