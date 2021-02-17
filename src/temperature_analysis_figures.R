############################################
# Models and plots for VOC Note
# TS: 2020-02-17
#

setwd("~/Documents/covid19-variant-N501Y/")

library(tidyverse)
library(gridExtra)
library(sf)
library(lme4)


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

summary(lm45)
summary(lm46)
summary(lm47)
summary(lm48)
summary(lm49)
summary(lm50)

# correct Rt by attack rate

lm45_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d45)
lm46_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d46)
lm47_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d47)
lm48_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d48)
lm49_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d49)
lm50_ar <- lm(Rtnext/(1-attackrate) ~ s.freq + (s.temperature + s.pop)*s.freq, d50)

summary(lm45_ar)
summary(lm46_ar)
summary(lm47_ar)
summary(lm48_ar)
summary(lm49_ar)
summary(lm50_ar)

# calculate differences in coefficients between models for each week
round(coefficients(lm45_ar) - coefficients(lm45), digits = 3)
round(coefficients(lm46_ar) - coefficients(lm46), digits = 3)
round(coefficients(lm47_ar) - coefficients(lm47), digits = 3)
round(coefficients(lm48_ar) - coefficients(lm48), digits = 3)
round(coefficients(lm49_ar) - coefficients(lm49), digits = 3)
round(coefficients(lm50_ar) - coefficients(lm50), digits = 3)

# get r-squared
round(-(summary(lm45)$r.squared - summary(lm45_ar)$r.squared), digits = 2)
round(-(summary(lm46)$r.squared - summary(lm46_ar)$r.squared), digits = 2)
round(-(summary(lm47)$r.squared - summary(lm47_ar)$r.squared), digits = 2)
round(-(summary(lm48)$r.squared - summary(lm48_ar)$r.squared), digits = 2)
round(-(summary(lm49)$r.squared - summary(lm49_ar)$r.squared), digits = 2)
round(-(summary(lm50)$r.squared - summary(lm50_ar)$r.squared), digits = 2)


# test effect of tiers
# (removing lone tier 1 isle of wight)
summary(lm(Rtnext ~ s.freq + as.factor(tier) + (s.temperature + s.pop)*s.freq, d49[d49$tier %in% c(2,3),]))
summary(lm(Rtnext ~ s.freq + as.factor(tier) + (s.temperature + s.pop)*s.freq, d50[d50$tier %in% c(2,3),]))

summary(lm(Rt ~ s.freq + as.factor(tier) + (s.temperature + s.pop)*s.freq, combined_data))

t.test(d50[d50$tier == 2,]$Rtnext, d50[d50$tier == 3,]$Rtnext)

# test ratios of R of VOC and non-VOC strains
# I feel like its valid keeping everything together here,
# because the ratio isn't something that should really vary temporally?
# i.e. the ratio shouldn't be effected by differences in NPIs between dates, etc.
transmission_lm <- lm(Ratio ~ as.factor(epiweek) + name + temperature, 
                      data = transmission_output_joint[transmission_output_joint$epiweek %in% c(45:50),])
summary(transmission_lm)
transmission_lm$coefficients["temperature"]
confint(transmission_lm, "temperature")

transmission_lmer <- lmer(Ratio ~ temperature + as.factor(epiweek) + (1|name), 
                         data=transmission_output_joint[transmission_output_joint$epiweek %in% c(45:50),])
lme4::fixef(transmission_lmer)["temperature"]
confint(transmission_lmer, "temperature")


####################################
# Figure Plotting

# make some extra fields for plotting legends
combined_data$week <- paste("Week", combined_data$epiweek)
combined_data$lockdown <- as.character(NA)
combined_data[combined_data$tier == 5,]$lockdown <- "National Lockdown"
combined_data[combined_data$tier %in% c(1,2,3,4),]$lockdown <- "Tier System"

# Big plot of national lockdown R vs Tiered lockdown R across the VOC spreading weeks
nat_lock <- ggplot(combined_data[combined_data$epiweek %in% c(43, 44, 45, 46, 47, 48, 49, 50),], 
                   aes(x = as.factor(epiweek), y = Rtnext)) + 
  geom_boxplot(outlier.size = 0, alpha = 0.5) + 
  geom_point(aes(fill = novel_frac*100), shape = 21, size = 3, position = position_jitterdodge(jitter.width = 0.4), alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Week",
       y = expression(R[t]),
       fill = "VOC %") +
  geom_vline(xintercept = 2.5) +
  geom_vline(xintercept = 6.5) +
  second_theme +
  theme(legend.position = c(0.32, 0.78))
nat_lock

# density plots of VOC distribution across those weeks
VOC_dist_full <- ggplot(combined_data[combined_data$epiweek %in% c(43, 44, 45, 46, 47, 48, 49, 50),],
                        aes(x = novel_frac*100)) +
  geom_density() +
  facet_wrap(~week, scales = "free_y", nrow = 1) +
  second_theme +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45))
VOC_dist_full

# effect of tiers in week 50
tier_plot <- ggplot(combined_data[combined_data$epiweek == 50 &
                                          combined_data$tier %in% c(2,3),], 
                          aes(x = as.factor(tier), y = Rtnext, fill = as.factor(tier))) + 
  geom_boxplot(outlier.size = 0, alpha = 0.5) +
  scale_fill_manual(values=c("grey20", "grey10", "grey0")) +
  geom_point(shape = 21, size = 3, position = position_jitterdodge(jitter.width = 0.8), alpha = 0.7) +
  labs(x = "Lockdown Tier",
       y = expression(R[t])) +
  second_theme +
  theme(legend.position = "none")
tier_plot


ggsave("figures/national_lockdown.svg", nat_lock, height = 5, width = 16)
ggsave("figures/VOC_dist_full.svg", VOC_dist_full, height = 2, width = 15)
ggsave("figures/Rt_tiers.svg", tier_plot, width = 4, height = 4)

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
  scale_fill_gradient(low = "blue", high = "yellow") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Temperature (°C)",
       y = "VOC Frequency",
       title = "Week 46, Beginning 2020-11-09",
       fill = expression(R[t])) +
  main_theme +
  theme(aspect.ratio = 1)
heatmap_plot_46


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
  scale_fill_gradient(low = "blue", high = "yellow") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Temperature (°C)",
       y = "VOC Frequency",
       title = "Week 50, Beginning 2020-12-07",
       fill = expression(R[t])) +
  main_theme +
  theme(aspect.ratio = 1)
heatmap_plot_50

# tiff because it seems to preserve colour contrast better than svg
ggsave("figures/week_46_independent_scale.tiff", heatmap_plot_46, height = 5, width = 5)
ggsave("figures/week_50_independent_scale.tiff", heatmap_plot_50, height = 5, width = 5)



#################################
# Supplementary plots

ratios_plot <- ggplot(transmission_output_joint[transmission_output_joint$epiweek %in% c(45:50),],
       aes(x = temperature, y = Ratio)) + geom_point(aes(col = as.factor(epiweek))) +
  geom_smooth(method = lm) +
  labs(x = "Temperature (°C)",
       y = expression(paste("Ratio VOC to non-VOC ", R[t])),
       col = "Week") +
  second_theme
ratios_plot

ggsave("figures/ratios_plot.svg", ratios_plot, height = 5, width = 7)
