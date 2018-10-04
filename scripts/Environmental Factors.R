# title: Summer 2017 Fungal Necromass Collaboration Environmental Factors
# author: KVB
#date: "9/21/2018"

# removing objects from the R Environment
rm(list=ls())

# load required packages
install.packages("pacman")
pacman::p_load(tidyverse,dplyr,broom, minpack.lm,Hmisc,car,gtools,lme4,nlme,minpack.lm)

# read in data file
fundecomp <- read.delim("data/CompiledNecroDecomp_summer2017.csv", sep = ",")

# subsetting data to just include unique env data for each site
envdata <- fundecomp  %>%
  select(site:harvest_date,-isolate) %>%
  filter(!is.na(soil_moisture)) %>%
  distinct()

# averages 
env.avg <- envdata %>%
  select(setting,treatment,pH,soil_moisture) %>%
  group_by(setting,treatment) %>%
  summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(env.avg, "output/avgs_envdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
