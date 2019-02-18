# title: Summer 2017 Fungal Necromass Collaboration Environmental Factors
# author: KVB
#date: "9/21/2018"

# removing objects from the R Environment
rm(list=ls())

# load required packages
packload <- c("tidyverse","dplyr","RColorBrewer", "broom","lmerTest","minpack.lm","Hmisc","car","gtools","lme4","nlme","minpack.lm","visreg", "lsmeans")
lapply(packload, library, character.only=TRUE)

# read in data file
fundecomp <- read.delim("data/CompiledNecroDecomp_summer2017.csv", sep = ",")

# subsetting data to just include unique env data for each site
envdata <- fundecomp  %>%
  select(site:myc_assoc,-plot_Id_pair, -unique_plotID,-unique_sampID,-treatment, incub_period.days) %>%
  filter(!is.na(soil_moisture)) %>%
  filter(myc_assoc != "mixed") %>%
  distinct()

envdata_mr <- fundecomp  %>%
  select(site,location,soil_moisture, pH, myc_assoc, isolate, mass_remaining) %>%
  filter(!is.na(mass_remaining)) %>%
  distinct()

envdata_mr.iso <- fundecomp  %>%
  select(site,location,soil_moisture, pH, myc_assoc, isolate, mass_remaining) %>%
  filter(!is.na(mass_remaining)) %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(myc_assoc != "mixed") %>%
  distinct()

MCf_env <- envdata   %>%
  filter(site !="Oak Savanna" & site !="Prairie") 

CCp_env <-  envdata   %>%
  filter(site !="Oak Savanna" & site !="Temperate Forest") 


CCos_env <-  envdata   %>%
  filter(site !="Prairie" & site !="Temperate Forest") 


# averages 
env.avg <- envdata %>%
  select(site,myc_assoc,pH,soil_moisture) %>%
  group_by(site,myc_assoc) %>%
  summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(env.avg, "output/avgs_envdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# environmental controls
sp_moist = lm(mass_remaining ~ soil_moisture, data=envdata_mr,na.action = na.exclude) 
summary(sp_moist)

sp_pH= lm(mass_remaining ~ pH, data=envdata_mr,na.action = na.exclude) 
summary(sp_moist)

# Plots
# seeing if pH/ soil moisture differ among mycorrhizal association/site
positions <- c("Prairie", "Oak Savanna", "Temperate Forest")
pH.site <- ggplot(envdata, aes(x=site, y=pH)) +
  geom_boxplot(width=0.5) +
  geom_jitter(size = 2.5, width = 0.2,alpha=0.5,aes(color = myc_assoc)) +
  scale_colour_manual(labels = c("AM", "EcM"), values = c("red", "blue")) +
  scale_x_discrete(limits = positions) +
  xlab("Site") + 
  ylab("pH") +
  labs(color = "Myc Association")  +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
ggsave(pH.site, filename = "./figures/pH.site.png", units = "in", width = 10 , height = 5, dpi = 300)


moist.site <- ggplot(envdata, aes(x=site, y=soil_moisture)) +
  geom_boxplot(width=0.5) +
  geom_jitter(size = 2.5, width = 0.2,alpha=0.5,aes(color = myc_assoc)) +
  scale_colour_manual(labels = c("AM", "EcM"), values = c("red", "blue")) +
  xlab("Site") + 
  ylab("Water Content") +
  labs(color = "Myc Association")  +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
ggsave(pH.site, filename = "./figures/pH.site.png", units = "in", width = 10 , height = 5, dpi = 300)

# Stats
# seeing if pH/ soil moisture differ among mycorrhizal association/site
# pH
pHmod <- lm(pH ~  myc_assoc * site, data=envdata, na.action=na.exclude)
summary(pHmod)
plot(pHmod) 
Anova(pHmod,type = "III")
qqnorm(resid(pHmod))
qqline(resid(pHmod))
hist(resid(pHmod))


# ls means comparisons vegetation type
ph.veg = lsmeans(pHmod, ~ site)
pairs(ph.veg,adjust="tukey")
phlph.vegpw = cld(ph.veg,alpha = 0.05,Letters = letters, adjust  = "tukey") 
phlph.vegpw 

# Plot
ggplot(phlph.vegpw, aes(x = site, y = lsmean, shape=myc_assoc, color = .group)) +
  geom_point(size  = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n pH")


# ls means comparisons vegetation type/ mycorrhizal association 
ph.mycveg = lsmeans(pHmod, ~ site + myc_assoc)
pairs(ph.mycveg,adjust="tukey")
phlph.mycvegpw = cld(ph.mycveg,alpha = 0.05,Letters = letters, adjust  = "tukey") 
phlph.mycvegpw 

# Plot
ggplot(phlph.mycvegpw, aes(x = site, y = lsmean, shape=myc_assoc, color = .group)) +
  geom_point(size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n pH")

# soil moisture 
mCmoist <- lm(log(soil_moisture) ~  myc_assoc * incub_period.days , data=MCf_env, na.action=na.exclude)
summary(mCmoist)
plot(mCmoist) 
Anova(mCmoist,type = "III")
qqnorm(resid(mCmoist))
qqline(resid(mCmoist))

# ls means comparisons vegetation type/ mycorrhizal association 
moist.mycveg = lsmeans(mCmoist, ~ myc_assoc)
pairs(moist.mycveg, adjust="tukey")
moist.mycvegpw  = cld(moist.mycveg,alpha = 0.05,Letters = letters, adjust  = "tukey") 
moist.mycvegpw

# Plot
ggplot(moist.mycvegpw ,aes(x = site, y = lsmean, shape=myc_assoc, color = .group)) +
  geom_point(size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax  =  upper.CL),width =  0.2, size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n moisture")

# soil moisture 
CCpmoist <- lm(log(soil_moisture) ~  myc_assoc  , data=CCp_env, na.action=na.exclude)
summary(CCpmoist)
plot(CCmoist) 
Anova(CCpmoist,type = "III")
qqnorm(resid(CCpmoist))
qqline(resid(CCmoist))

# soil moisture 
CCosmoist <- lm(log(soil_moisture) ~  myc_assoc * incub_period.days , data=CCos_env, na.action=na.exclude)
summary(CCosmoist)
plot(CCosmoist) 
Anova(CCosmoist,type = "III")
qqnorm(resid(CCosmoist))
qqline(resid(CCosmoist))