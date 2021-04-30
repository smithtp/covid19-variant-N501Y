############################################
# Combine Rt, VOC frequency and Ratios data
# with climate and population variables
# TS: 2020-02-17
#

library(tidyverse)
library(sf)

### prepLTLA 
la_rt<- readRDS("data/la_rt.rds")%>%
  group_by(area)%>%mutate(Rtnext=lead(Rt,order_by=epiweek),
                          CIupnext=lead(CIup,order_by=epiweek),
                          CIlownext=lead(CIlow,order_by=epiweek)
  )%>%ungroup

LTLA_sgss <- readRDS("data/sgss_la.rds")%>%
  mutate(novel_frac=sgss_s_negative_corrected/(sgss_s_negative_corrected+sgss_s_positive_corrected))%>%
  inner_join(la_rt,by=c("area","epiweek"))%>%
  dplyr::rename(novel=sgss_s_negative,other=sgss_s_positive)

#############################
# ---  Add Temperature  --- #
#############################

uk_temperature_long <- read.csv("data/UK-LTLA-temperature.csv")
uk_temp <- uk_temperature_long[,c("area", "temperature", "week")] %>%
  group_by(area = area, epiweek = week) %>%
  dplyr::summarize(temperature = mean(temperature, na.rm=TRUE))

# merge temperature data with Rt data
combined_data <- left_join(LTLA_sgss, uk_temp, by=c("epiweek","area"))

#############################
# ---  Add Attack Rate  --- #
#############################

# case data to calculate attack rates
# note, these are estimated infections from the model - not real infections?
ltla_rt_cases <- readRDS("data/aggregates_infections_rt.rds") # swapnil aggregated these with the cases
ltla_rt_cases$epiweek <- as.numeric(format(as.Date(ltla_rt_cases$period_start), "%V"))
ltla_rt_cases$epiyear <- as.numeric(format(as.Date(ltla_rt_cases$period_start), "%Y"))

ltla_cases <- ltla_rt_cases[ltla_rt_cases$type == "Infections" &
                              ltla_rt_cases$epiyear == 2020,] %>%
  group_by(area, epiweek) %>%
  dplyr::summarize(infections = sum(value)) %>%
  dplyr::mutate(cumulativeinfections = cumsum(infections))

# add total population to calculate attack rate
population_data <- read.csv("data/modified_population.csv") %>%
  dplyr::rename(area = AREA, population = Y2018)
ltla_cases <- left_join(ltla_cases, population_data[,c("area", "population")], by = "area") %>%
  mutate(attackrate = cumulativeinfections/population, 
         infections_per_100k = infections*(100000/population))


combined_data <- left_join(combined_data, ltla_cases, by=c("epiweek","area"))

####################################
# ---  Add Population density  --- #
####################################


uk_popdensity <- data.frame(readRDS("data/population-density-UK-LTLA.RDS")) 
uk_popdensity$area <- rownames(uk_popdensity)

combined_data <- left_join(combined_data, uk_popdensity, by = "area")

################################
# ---  Add Lockdown Tiers  --- #
################################


# add the UK tiers
uk_tiers_data <- read.csv("data/NPI_created_dataset_20_01_2021.csv")

# shorten to just the columns we want
uk_tiers_data <- uk_tiers_data[,c("date", "ltla", "tier_1", "tier_2", "tier_3", "tier_4", "national_lockdown")]

# turn tiers into single column
uk_tiers <- uk_tiers_data %>%
  group_by(ltla = ltla, date = date) %>%
  gather(tier, val, -date, -ltla) %>%
  filter(val == 1) %>%
  dplyr::select(-val) %>%
  dplyr::arrange(ltla, date)

# merge this back into the original tier data to get NAs where there was no tier
uk_tiers_data <- left_join(uk_tiers_data, uk_tiers, by = c("date", "ltla"))

# and change to numeric tiers for later regressions
uk_tiers_data[is.na(uk_tiers_data$tier),]$tier <- 0
uk_tiers_data[uk_tiers_data$tier == "tier_1",]$tier <- 1
uk_tiers_data[uk_tiers_data$tier == "tier_2",]$tier <- 2
uk_tiers_data[uk_tiers_data$tier == "tier_3",]$tier <- 3
uk_tiers_data[uk_tiers_data$tier == "tier_4",]$tier <- 4
uk_tiers_data[uk_tiers_data$tier == "national_lockdown",]$tier <- 5
uk_tiers_data$tier <- as.numeric(uk_tiers_data$tier)

uk_tiers_data$epiweek <- format(as.Date(uk_tiers_data$date, format='%d/%m/%Y'), "%V")
uk_tiers_data$epiyear <- format(as.Date(uk_tiers_data$date, format='%d/%m/%Y'), "%Y")


# group by week
# need a mode function
mode <- function(codes){
  as.numeric(which.max(tabulate(codes)))
}

weekly_tiers <- uk_tiers_data[,c("ltla", "tier", "epiweek", "epiyear")] %>%
  group_by(area = ltla, epiweek = epiweek, epiyear = epiyear) %>%
  dplyr::summarize(tier = mode(tier))

# drop to just 2020 for now
weekly_tiers <- weekly_tiers[weekly_tiers$epiyear == "2020",]
weekly_tiers$epiweek <- as.numeric(weekly_tiers$epiweek)
weekly_tiers$area <- as.character(weekly_tiers$area)

# merge into Rt data
combined_data <- left_join(combined_data, weekly_tiers[,c("area", "epiweek", "tier")], by = c("area", "epiweek"))


### --- DONE!

#####################
# seperately, we'd like to compare the ratios of Rts
# of VOC and non-VOC strains with temperature -
# are there temperatures where VOC does better than non-VOC for example?
# but, these ratios are at NHS region level, rather than LTLA
# "NHS England Sustainability and Transformation Plan (STP) areas"

#transmission_output_joint <- readRDS("data/transmission_output_joint.rds")

# updated data from Swapnil:
STP_rts <- read.csv("data/df_all_orig_tg.csv")
names(STP_rts) <- c("name", "epiweek", "R(S-)", "R(S+)")

STP_rts$Ratio <- STP_rts$`R(S-)`/STP_rts$`R(S+)`

# get epiweek into sensible format with gsub
# transmission_output_joint$epiweek <- as.numeric(gsub('\\D+','', transmission_output_joint$epiweek))
STP_rts$epiweek <- as.numeric(gsub('\\D+','', STP_rts$epiweek))

# never fear, I've used the geometry here to average across climate regions
stp_temperature <- as.data.frame(readRDS("data/temp-UK-STP.RDS"))

make_long <- function(df, clim_var){
  df$area <- row.names(df)
  return(pivot_longer(df,
                      cols = c(1:(ncol(df)-1)),
                      names_to = "date",
                      values_to = clim_var))
}

stp_temperature_long <- make_long(stp_temperature, "temperature")
stp_temperature_long$week <- as.numeric(format(as.Date(stp_temperature_long$date), "%V"))

stp_temp <- stp_temperature_long[,c("area", "temperature", "week")] %>%
  group_by(name = area, epiweek = week) %>%
  dplyr::summarize(temperature = mean(temperature, na.rm=TRUE))

# merge with transmission data
# transmission_output_joint <- left_join(transmission_output_joint, stp_temp, by =  c("name", "epiweek"))
STP_rts <- left_join(STP_rts, stp_temp, by =  c("name", "epiweek"))


### --- DONE!

# save the outputs
saveRDS(combined_data, "data/combine_rt_climate.rds")
saveRDS(STP_rts, "data/combine_transmission_output_climate.rds")
