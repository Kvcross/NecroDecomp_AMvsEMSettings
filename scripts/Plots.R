# title: Summer 2017 Fungal Necromass Collaboration Plots
# author: KVB
#date: "9/25/2018"

# housekeeping
rm(list=ls())

# load required packages
install.packages("pacman")
pacman::p_load(tidyverse,dplyr)

# read in data file- with intitial values for MR
ifundecomp <- read.delim("data/NecroDecompinitial_summer2017.csv", sep = ",")

# removing NAs from mass remaining
ifundecomp <- ifundecomp %>% 
  filter(!is.na(mass_remaining))

# subsetting
ifundecomp_myc <- ifundecomp %>%
  filter(!myc_assoc %in% "mixed")

ifundecomp_melb.mort <- ifundecomp_myc %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(!myc_assoc %in% "mixed")

# mass Remaining for each isolate for the different incubation periods
MR1 <-  ggplot(ifundecomp, aes(x=incub_period.days, y=mass_remaining, color=isolate)) +
  geom_point() +
  scale_color_manual(values=c("black", "dimgrey", "gray","cornsilk3")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(setting ~.)

MR2 <- ggplot(ifundecomp, aes(x=incub_period.days, y=mass_remaining, color=mel_status)) +
  geom_point() +
  scale_color_manual(values=c("black", "dimgrey", "gray","cornsilk3")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(setting ~.)


MR3 <-ggplot(ifundecomp, aes(x=incub_period.days, y=mass_remaining, color=isolate)) +
  stat_summary(fun.data = "mean_se", size = 0.5) +
  stat_summary(fun.y = mean,geom = "line", aes(group=isolate), size=0.5) +
  stat_summary(fun.y="mean", geom="point", size=0.5, 
               aes(group = isolate), show.legend=FALSE) +
  scale_color_manual(values=c("black", "dimgrey", "gray", "cornsilk3", "lemonchiffon4")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(.~setting,scales = "free")
  
# mass Remaining for each isolate for the different incubation periods by site
MR4 <-  ggplot(ifundecomp_melb.mort, aes(x=incub_period.days, y=mass_remaining, color=isolate)) +
  stat_summary(fun.data = "mean_se", size = 0.5) +
  stat_summary(fun.y = mean,geom = "line", aes(group=isolate), size=0.5) +
  stat_summary(fun.y="mean", geom="point", size=0.5, 
               aes(group = isolate), show.legend=FALSE) +
  scale_color_manual(values=c("black", "dimgrey")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(.~setting,scales = "free")

# mass Remaining for each isolate for the different incubation periods by site/ myc assoc
MR5 <-  ggplot(ifundecomp_melb.mort, aes(x=incub_period.days, y=mass_remaining, color=isolate, shape=isolate_myc)) +
  stat_summary(fun.data = "mean_se", size = 0.5) +
  stat_summary(fun.y = mean,geom = "line", aes(group=isolate_myc), size=0.5) +
  stat_summary(fun.y="mean", geom="point", size=0.5, 
               aes(group = isolate_myc), show.legend=FALSE) +
  scale_color_manual(values=c("black", "dimgrey")) +
  scale_shape_manual(name = "Treatment & Isolate", values = c( 1,1,19,19)) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(.~setting,scales = "free")

# pH plot
pH_set = ggplot(ifundecomp, aes(x=pH, y=mass_remaining)) + 
  geom_point(color="grey") +
  geom_smooth(method=lm, se=FALSE, size=0.5, color="black") +
  theme_classic(base_size = 20) +
  ylim(0,50) +
  ylab("Mass Remaining (%)") +
  xlab("pH") +
  labs(color = "Fungal Species") 

# moisture plot
moist_set = ggplot(ifundecomp, aes(x=soil_moisture, y=mass_remaining)) + 
  geom_point(color="grey") +
  geom_smooth(method=lm, se=FALSE, size=0.5, color="black") +
  theme_classic(base_size = 20) +
  ylim(0,50) +
  ylab("Mass Remaining (%)") +
  xlab("Volumetric Water Content (%)") +
  labs(color = "Fungal Species") 


ggsave(MR3, filename = "./figures/Prelim_MassRemain_sp.png", units = "in", width = 15, height = 6, dpi = 300)
ggsave(MR4, filename = "./figures/Prelim_MassRemain_mel.png", units = "in", width = 15 , height = 6, dpi = 300)
ggsave(MR5, filename = "./figures/Prelim_MassRemain_myc.png", units = "in", width = 15 , height = 6, dpi = 300)
ggsave(pH_set, filename = "./figures/MassRemain_pH.png", units = "in", width = 9, height = 7, dpi = 300)
ggsave(moist_set, filename = "./figures/MassRemain_moist.png", units = "in", width = 9, height = 7, dpi = 300)
ggsave(pH_set, filename = "./figures/MassRemain_pH.png", units = "in", width = 9, height = 7, dpi = 300)

