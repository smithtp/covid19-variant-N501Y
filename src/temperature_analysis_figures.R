############################################
# Models and plots for VOC Note
# TS: 2020-02-17
#

library(tidyverse)
library(lme4)
library(xtable)

######################
# Plotting themes:

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size=16, vjust=1),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

second_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size=20, vjust=1),
        legend.text=element_text(size=20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 18))

# function to calculate Rt given a particular model of VOC/temp/pop
calc_Rt <- function(model, VOC_freq, temperature, log_pop_density){
  # pull out model coefficients
  intercept <- model$coefficients[[1]]
  b1 <- model$coefficients[[2]]
  b2 <- model$coefficients[[3]]
  b3 <- model$coefficients[[4]]
  b4 <- model$coefficients[[5]]
  b5 <- model$coefficients[[6]]
  # equation:
  Rt <- intercept + b1*VOC_freq + b2*temperature + b3*log_pop_density + b4*VOC_freq*temperature +
    b5*VOC_freq*log_pop_density
  return(Rt)
}

#############################

# read data

combined_data <- readRDS("data/combine_rt_climate.rds")
transmission_output_joint <- readRDS("data/combine_transmission_output_climate.rds")

# scale variables
combined_data$s.freq <- scale(combined_data$novel_frac)
combined_data$s.temperature <- scale(combined_data$temperature)
combined_data$log.pop <- log10(combined_data$Pop_density)
combined_data$s.pop <- scale(log10(combined_data$Pop_density))


#######################################
# Analysis

# pull out each week and test for the effect of the new variant and environment
d45 <- combined_data[combined_data$epiweek == 45,]
d46 <- combined_data[combined_data$epiweek == 46,]
d47 <- combined_data[combined_data$epiweek == 47,]
d48 <- combined_data[combined_data$epiweek == 48,]
d49 <- combined_data[combined_data$epiweek == 49,]
d50 <- combined_data[combined_data$epiweek == 50,]

lm45 <- lm(Rtnext ~ s.freq + (s.temperature + s.pop)*s.freq, d45)
lm46 <- lm(Rtnext ~ s.freq + (s.temperature + s.pop)*s.freq, d46)
lm47 <- lm(Rtnext ~ s.freq + (s.temperature + s.pop)*s.freq, d47)
lm48 <- lm(Rtnext ~ s.freq + (s.temperature + s.pop)*s.freq, d48)
lm49 <- lm(Rtnext ~ s.freq + (s.temperature + s.pop)*s.freq, d49)
lm50 <- lm(Rtnext ~ s.freq + (s.temperature + s.pop)*s.freq, d50)

table_1 <- function(){
  print("Table 1: VOC, Temperature, Pop density vs R with interactions week 46 regression model with scaled coefficients:")
  return(xtable(summary(lm46)))
}

table_2 <- function(){
  print("Table 2: VOC, Temperature, Pop density vs R with interactions week 50 regression model with scaled coefficients:")
  return(xtable(summary(lm50)))
}

table_S1 <- function(){
  week49 <- summary(lm(Rtnext ~ s.freq + as.factor(tier) + (s.temperature + s.pop)*s.freq, d49[d49$tier %in% c(2,3),]))
  week50 <- summary(lm(Rtnext ~ s.freq + as.factor(tier) + (s.temperature + s.pop)*s.freq, d50[d50$tier %in% c(2,3),]))
  print("Table S1: VOC, Tier, Temperature, Pop density vs R with interactions week 49 and 50 regression models:")
  return(list(week49, week50))
}


table_S2 <- function(){
  # correct Rt by attack rate
  
  lm45_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d45)
  lm46_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d46)
  lm47_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d47)
  lm48_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d48)
  lm49_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d49)
  lm50_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d50)
  
  # calculate differences in coefficients between models for each week
  ar_df <- data.frame(rbind(round(coefficients(lm45_ar) - coefficients(lm45), digits = 3),
                            round(coefficients(lm46_ar) - coefficients(lm46), digits = 3),
                            round(coefficients(lm47_ar) - coefficients(lm47), digits = 3),
                            round(coefficients(lm48_ar) - coefficients(lm48), digits = 3),
                            round(coefficients(lm49_ar) - coefficients(lm49), digits = 3),
                            round(coefficients(lm50_ar) - coefficients(lm50), digits = 3)))
  names(ar_df) <- c("Delta Intercept", "Delta VOC", "Delta Temp", "Delta Pop", "Delta VOCxTemp", "Delta VOCxPop")

  ar_df$Week <- seq(45, 50, 1)

  # get r-squared
  ar_df$`Delta r2` <- c(round(-(summary(lm45)$r.squared - summary(lm45_ar)$r.squared), digits = 2),
                        round(-(summary(lm46)$r.squared - summary(lm46_ar)$r.squared), digits = 2),
                        round(-(summary(lm47)$r.squared - summary(lm47_ar)$r.squared), digits = 2),
                        round(-(summary(lm48)$r.squared - summary(lm48_ar)$r.squared), digits = 2),
                        round(-(summary(lm49)$r.squared - summary(lm49_ar)$r.squared), digits = 2),
                        round(-(summary(lm50)$r.squared - summary(lm50_ar)$r.squared), digits = 2))
  print("Table S2: Difference in model coefficients and r-squared when using attack rate-corrected Rt valules:")
  return(xtable(ar_df))
}

# test ratios of R of VOC and non-VOC strains
# I feel like its valid keeping everything together here,
# because the ratio isn't something that should really vary temporally?
# i.e. the ratio shouldn't be effected by differences in NPIs between dates, etc.
table_S3 <- function(){
  
  transmission_lm <- lm(Ratio ~ as.factor(epiweek) + name + temperature, 
                        data = transmission_output_joint[transmission_output_joint$epiweek %in% c(45:50),])
  
  transmission_lmer <- lmer(Ratio ~ temperature + as.factor(epiweek) + (1|name), 
                            data=transmission_output_joint[transmission_output_joint$epiweek %in% c(45:50),])
  
  model_df <- data.frame(Model = c("Fixed", "Random"),
                         Temeperature_coef = c(transmission_lm$coefficients["temperature"],
                                               lme4::fixef(transmission_lmer)["temperature"]),
                         Temp_lower_CI = c(confint(transmission_lm, "temperature")[1],
                                           confint(transmission_lmer, "temperature")[1]),
                         Temp_upper_CI = c(confint(transmission_lm, "temperature")[2],
                                           confint(transmission_lmer, "temperature")[2]))
  
  print("Table S3: Effects of temperature on Ratio of Rt - fixed and random effects models:")
  return(xtable(model_df))
}

####################################
# Figure Plotting

# make some extra fields for plotting legends
combined_data$week <- paste("Week", combined_data$epiweek)

combined_data <- combined_data %>%
  group_by(epiweek) %>%
  mutate(med_VOC = median(novel_frac))

# Big plot of national lockdown R vs Tiered lockdown R across the VOC spreading weeks
nat_lock <- ggplot(combined_data[combined_data$epiweek %in% c(43, 44, 45, 46, 47, 48, 49, 50),], 
                   aes(x = as.factor(epiweek), y = Rtnext)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = med_VOC*100)) + 
  geom_point(aes(fill = novel_frac*100), shape = 21, size = 3, position = position_jitterdodge(jitter.width = 0.4)) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Week",
       y = expression(R[t]),
       fill = "VOC %") +
  geom_vline(xintercept = 2.5) +
  geom_vline(xintercept = 6.5) +
  second_theme +
  theme(legend.position = c(0.32, 0.78)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())


# density plots of VOC distribution across those weeks
VOC_dist_full <- ggplot(combined_data[combined_data$epiweek %in% c(43, 44, 45, 46, 47, 48, 49, 50),],
                        aes(x = novel_frac*100)) +
  geom_density() +
  facet_wrap(~week, scales = "free_y", nrow = 1) +
  second_theme +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45))


# heatmap plots of model results in early and late VOC spread
unscaled_week46 <- lm(Rtnext ~ novel_frac + (temperature + log.pop)*novel_frac, d46)

temps <- seq(9.5, 13, by = 0.1)
freqs <- 0.01*(seq(-3, 100, by = 0.1))
pops <- mean(d46$log.pop) # just use median pop density?
grid <- expand.grid(temps, freqs, pops)
predicted_Rt_46 <- setNames(data.frame(grid), c("Temperature", "Frequency", "Pop_density"))

predicted_Rt_46$Rt <- calc_Rt(unscaled_week46, VOC_freq = predicted_Rt_46$Frequency,
                              temperature = predicted_Rt_46$Temperature,
                              log_pop_density = predicted_Rt_46$Pop_density)


heatmap_plot_46 <- ggplot(predicted_Rt_46, aes(x = Temperature, y = Frequency)) + 
  geom_tile(aes(fill = Rt)) +
  geom_point(data = d46, aes(x = temperature, y = novel_frac, fill = Rtnext), size = 4, shape = 21) +
  #scale_fill_viridis_c(limits = c(0.5, 3)) +
  scale_fill_gradient(low = "blue", high = "yellow", limits = c(0.55, 2.8)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Temperature (°C)",
       y = "VOC Frequency",
       title = "Week 46, Beginning 2020-11-09",
       fill = expression(R[t])) +
  main_theme +
  theme(aspect.ratio = 1,
        legend.position = "none")

# # with VOC frequency as the colour?
# ggplot(predicted_Rt_46, aes(x = Temperature, y = Rt)) + 
#   #geom_tile(aes(fill = Frequency)) +
#   geom_point(data = d46, aes(x = temperature, y = Rtnext, fill = novel_frac), size = 4, shape = 21) +
#   #scale_fill_viridis_c(limits = c(0.5, 3)) +
#   scale_fill_gradient(low = "blue", high = "yellow") +
#   #scale_x_continuous(expand = c(0, 0)) +
#   #scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "Temperature (°C)",
#        y = expression(R[t]),
#        title = "Week 46, Beginning 2020-11-09",
#        fill = "VOC Frequency") +
#   main_theme +
#   theme(aspect.ratio = 1)


unscaled_week50 <- lm(Rtnext ~ novel_frac + (temperature + log.pop)*novel_frac, d50)

temps <- seq(5, 10, by = 0.1)
freqs <- 0.01*(seq(-3, 105, by = 0.1))
pops <- mean(d50$log.pop) # just use median pop density?
grid <- expand.grid(temps, freqs, pops)

predicted_Rt_50 <- setNames(data.frame(grid), c("Temperature", "Frequency", "Pop_density"))

predicted_Rt_50$Rt <- calc_Rt(unscaled_week50, VOC_freq = predicted_Rt_50$Frequency, 
                              temperature = predicted_Rt_50$Temperature, 
                              log_pop_density = predicted_Rt_50$Pop_density)


heatmap_plot_50 <- ggplot(predicted_Rt_50, aes(x = Temperature, y = Frequency)) + 
  geom_tile(aes(fill = Rt)) +
  geom_point(data = d50, aes(x = temperature, y = novel_frac, fill = Rtnext), size = 4, shape = 21) +
  #scale_fill_viridis_c(limits = c(0.5, 3)) +
  scale_fill_gradient(low = "blue", high = "yellow", limits = c(0.55, 2.8)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Temperature (°C)",
       y = "VOC Frequency",
       title = "Week 50, Beginning 2020-12-07",
       fill = expression(R[t])) +
  main_theme +
  theme(aspect.ratio = 1,
        legend.position = "none")


# # with VOC frequency as the colour?
# ggplot(predicted_Rt_50, aes(x = Temperature, y = Rt)) + 
#   #geom_tile(aes(fill = Frequency)) +
#   geom_point(data = d50, aes(x = temperature, y = Rtnext, fill = novel_frac), size = 4, shape = 21) +
#   #scale_fill_viridis_c(limits = c(0.5, 3)) +
#   scale_fill_gradient(low = "blue", high = "yellow") +
#   #scale_x_continuous(expand = c(0, 0)) +
#   #scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "Temperature (°C)",
#        y = expression(R[t]),
#        title = "Week 50, Beginning 2020-12-07",
#        fill = "VOC Frequency") +
#   main_theme +
#   theme(aspect.ratio = 1)


#################################
# Supplementary plots

attackrates_plot <- ggplot(combined_data[combined_data$epiweek %in% c(45, 46, 47, 48, 49, 50),]) + 
  geom_density(aes(x=attackrate*100, colour=week),show.legend=FALSE, size = 2) +
  stat_density(aes(x=attackrate*100, colour=week),
               geom="line",position="identity", size = 0) + 
  guides(colour = guide_legend(override.aes=list(size=1))) +
  scale_colour_manual(values = c("black", "grey45", "darkslateblue", "cyan4", "orange2", "coral3")) +
  geom_vline(xintercept = median(d45$attackrate*100), col = "black", linetype = "dashed", lwd = 1) +
  geom_vline(xintercept = median(d50$attackrate*100), col = "coral3", linetype = "dashed", lwd = 1) +
  labs(x = "Attack Rate (%)",
       y = "Density") +
  second_theme +
  theme(legend.position = c(0.8, 0.6),
        legend.title = element_blank())


#######################
####      MAIN     ####
#######################


# Run analyses and generate latex tables
print("Week 46 analysis:")
table_1()
print("")
print("=============================================================")
print("Week 50 analysis:")
table_2()
print("")
print("=============================================================")
print("Effect of lockdown tiers (supplementary):")
table_S1()
print("")
print("=============================================================")
print("Attack Rates sensitivity analysis (supplementary):")
table_S2()
print("")
print("=============================================================")
print("Effects of temperature on ratio of VOC to non-VOC Rt (supplementary):")
table_S3()
print("")
print("=============================================================")


ggsave("figures/national_lockdown.svg", nat_lock, height = 5, width = 16)
ggsave("figures/VOC_dist_full.svg", VOC_dist_full, height = 2, width = 15)

# tiff because it seems to preserve colour contrast better than svg
ggsave("figures/week_46_same_scale.tiff", heatmap_plot_46, height = 5, width = 5)
ggsave("figures/week_50_same_scale.tiff", heatmap_plot_50, height = 5, width = 5)

ggsave("figures/attackrates_plot.tiff", attackrates_plot, height = 5, width = 7, dpi=300, compression = "lzw")
