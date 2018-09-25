# title: Summer 2017 Fungal Necromass Mass Remaining Collaboration
# author: KVB
 #date: "9/21/2018"

# removing objects from the R Environment
rm(list=ls())

# load required packages
install.packages("pacman")
pacman::p_load(tidyverse,dplyr,broom, minpack.lm,Hmisc,car,gtools,nlme)

# read in data file
MR_data <- read.delim("data/CompiledNecroDecomp_summer2017.csv", sep = ",")

# specifying factors for data file
MR_data <-  within(MR_data, {
            site = factor(site, ordered=FALSE)
            setting = factor(setting, ordered=FALSE)
            exp = factor(exp, ordered=FALSE)
            location = factor(location, ordered=FALSE)
            plot_Id = factor(plot_Id, ordered=FALSE)
            treatment = factor(treatment, ordered=FALSE)
            myc_assoc  = factor(myc_assoc, ordered=FALSE)
            isolate = factor(isolate, ordered=FALSE)
})

ENV_data <- MR_data %>%
            select(site:harvest_date,-isolate) %>%
            filter(!is.na(soil_moisture)) %>%
            distinct()

# averages 
env.avg <- ENV_data %>%
  select(site,setting,treatment,pH,soil_moisture) %>%
  group_by(site,setting,treatment) %>%
  summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(env.avg, "output/avgs_envdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

mass.avg = MR_data %>%
  select(incub_period.days,isolate, mass_remaining) %>%
  group_by(incub_period.days,isolate) %>%
  summarise_all(funs(length,mean(., na.rm = TRUE),sd(., na.rm = TRUE),se=sd(., na.rm = TRUE)/sqrt(n())))
write.table(mass.avg , "output/avgs_massdata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# preliminary plots

# mass Remaining for each isolate for the different incubation periods
MR1 = ggplot(MR_data, aes(x=incub_period.days, y=mass_remaining, color=isolate)) +
  stat_summary(fun.data = "mean_se", size = 0.5) +
  stat_summary(fun.y = mean,geom = "line", aes(group=isolate), size=0.5) +
  stat_summary(fun.y="mean", geom="point", size=0.5, 
               aes(group = isolate), show.legend=FALSE) +
  scale_color_manual(values=c("black", "dimgrey", "gray", "cornsilk3", "lemonchiffon4")) +
  scale_x_continuous("Time (days)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) 

# mass Remaining for each isolate for the different incubation periods by site
MR2 = ggplot(MR_data, aes(x=incub_period.months, y=mass_remaining, color=isolate)) +
  stat_summary(fun.data = "mean_se", size = 0.5) +
  stat_summary(fun.y = mean,geom = "line", aes(group=isolate), size=0.5) +
  stat_summary(fun.y="mean", geom="point", size=0.5, 
               aes(group = isolate), show.legend=FALSE) +
  scale_color_manual(values=c("black", "dimgrey", "gray", "cornsilk3", "lemonchiffon4")) +
  scale_x_continuous("Time (months)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_wrap(~setting,scales = "free")

# same plot at MR2 with myc assoc
MR3 = ggplot(MR_data, aes(x=incub_period.days, y=mass_remaining, color=isolate, shape=myc_assoc)) +
  stat_summary(fun.data = "mean_se", size = 0.5) +
  stat_summary(fun.y = mean,geom = "line", aes(group=isolate), size=0.5) +
  stat_summary(fun.y="mean", geom="point", size=0.5, 
               aes(group = isolate), show.legend=FALSE) +
  scale_color_manual(values=c("black", "dimgrey", "gray", "cornsilk3", "lemonchiffon4")) +
  scale_x_continuous("Time (months)") + 
  scale_y_continuous(name="Mass Remaining (%)", limits=c(0, 100)) +
  labs(color="Fungal Species") +
  theme_classic(base_size=20) +
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
  facet_wrap(~setting,scales = "free")

ggsave(MR1, filename = "./figures/Prelim_MassRemain_sp.png", units = "in", width = 12, height = 5, dpi = 300)
ggsave(MR2, filename = "./figures/Prelim_MassRemain_site.png", units = "in", width = 12, height = 5, dpi = 300)
ggsave(MR3, filename = "./figures/MassRemain_site_myc.png", units = "in", width = 12, height = 5, dpi = 300)
