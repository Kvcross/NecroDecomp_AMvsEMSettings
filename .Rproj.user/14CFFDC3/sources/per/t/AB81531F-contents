# title: Summer 2017 Fungal Necromass Collaboration Plots
# author: KVB
#date: "9/25/2018"

# housekeeping
rm(list=ls())

# load required packages
packload <- c("tidyverse","dplyr","broom", "ggplot2", "Rmisc")
lapply(packload, library, character.only=TRUE)

# read in data file- with intitial values for MR
ifundecomp <- read.delim("data/NecroDecompinitial_summer2017.csv", sep = ",")

# removing NAs from mass remaining
ifundecomp <- ifundecomp %>% 
  filter(!is.na(mass_remaining))

# subsetting
ifundecomp_myc <- ifundecomp %>%
  filter(!myc_assoc %in% "mixed") 
  
ifundecomp_mb.me <- ifundecomp_myc %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(!myc_assoc %in% "mixed")

ifundecomp_mb.me.f5 <- ifundecomp_myc %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(!myc_assoc %in% "mixed") %>% 
  filter(!location %in% c("field_0","field_80"))

ifundecomp_mb.me.p2 <- ifundecomp_myc %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(!myc_assoc %in% "mixed") %>% 
  filter(!exp %in% "CC_praire1")

ifundecomp_mb.me.p2.early <- ifundecomp_mb.me.p2 %>%
  filter(!isolate %in% c("Ceno","Mel_White")) %>%
  filter(!myc_assoc %in% "mixed") %>% 
  filter(!exp %in% "CC_praire1")  %>% 
  filter(incub_period.days %in% c("11","14"))

# mass Remaining for each isolate for the different incubation periods
MR1 <-  ggplot(ifundecomp, aes(x=incub_period.days, y=mass_remaining, color=isolate)) +
  geom_point() +
  geom_smooth() +
  scale_color_manual(values=c("black", "dimgrey", "gray","cornsilk3")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(site ~.)

MR2 <- ggplot(ifundecomp, aes(x=incub_period.days, y=mass_remaining, color=mel_status)) +
  geom_point() +
  scale_color_manual(values=c("black", "dimgrey", "gray","cornsilk3")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(site ~.)
  

# mass Remaining for each isolate for the different incubation periods by site/ myc assoc
ifundecomp_mb.me$facet = factor(ifundecomp_mb.me$site, levels = c("Prairie", "Oak Savanna", "Temperate Forest"))
MR3 <-  ggplot(ifundecomp_mb.me, aes(x=incub_period.wks, y=mass_remaining,group=isolate_myc, shape=isolate_myc)) +
         stat_summary(fun.data = "mean_se", geom = "errorbar",size=1.1, color="gray25", show.legend=FALSE) +
         stat_summary(fun.y="mean", geom="line", size=1,color="gray35", show.legend=FALSE) +
         stat_summary(fun.data = "mean_se", geom="point", size=8,  stroke = 1.2, color="gray25",aes(shape = isolate_myc), show.legend=FALSE) +
         scale_shape_manual(name = "Treatment & Isolate",values = c(19,1,17,2)) +
         scale_x_continuous("Time (weeks)") + 
         scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
         theme_classic(base_size=30) +
         theme(legend.key=element_blank()) +
         theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
         facet_grid(.~facet,scales = "free")

ggsave(MR3, filename = "./figures/Prelim_MassRemain_myc.png", units = "in", width = 20 , height = 9, dpi = 300)

# Leaving out prairie exp 1: mass Remaining for each isolate for the different incubation periods by site/ myc assoc
ifundecomp_mb.me.p2$facet = factor(ifundecomp_mb.me.p2$site, levels = c("Prairie", "Oak Savanna", "Temperate Forest"))

MR4 <-  ggplot(ifundecomp_mb.me.p2, aes(x=incub_period.wks, y=mass_remaining,group=isolate_myc, shape=isolate_myc)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",size=1.1, color="gray25", show.legend=FALSE) +
  stat_summary(fun.y="mean", geom="line", size=1,color="gray35", show.legend=FALSE) +
  stat_summary(fun.data = "mean_se", geom="point", size=8,  stroke = 1.2, color="gray25",aes(shape = isolate_myc), show.legend=FALSE) +
  scale_shape_manual(name = "Treatment & Isolate",values = c(19,1,17,2)) +
  scale_x_continuous("Time (weeks)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  theme_classic(base_size=30) +
  theme(legend.key=element_blank()) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_grid(.~facet,scales = "free")
ggsave(MR4, filename = "./figures/Prelim_MassRemain_mycPr2.png", units = "in", width = 20 , height = 9, dpi = 300)

# Early Decay Stage 
MRE <- summarySE(ifundecomp_mb.me.p2.early, measurevar="mass_remaining", groupvars=c("site","isolate"))
positions <- c("Prairie", "Oak Savanna", "Temperate Forest")

MR5 <- ggplot(MRE, aes(x=site, y=mass_remaining, fill=isolate)) + 
       geom_bar(position=position_dodge(), stat="identity", color="black") +
       geom_errorbar(aes(ymin=mass_remaining-se, ymax=mass_remaining+se),width=.2,                   
                position=position_dodge(.9)) +
      scale_fill_manual(name="", labels = c("M. bicolor", "M. elongata"), values=c("gray35", "white"))+
       scale_x_discrete(limits = positions) +
       ylim(0,65) +
       ylab("Mass Remaining (%)") +
       xlab("Site") +
       labs(color = "Fungal Species") +
       theme_classic(base_size = 25)
ggsave(MR5, filename = "./figures/MRearly.png", units = "in", width = 11 , height = 9, dpi = 300)

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



