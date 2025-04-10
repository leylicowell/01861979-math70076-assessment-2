# in this script, we implement a bayesian regression model for economical damages
# using our merged data set

library(data.table) 
library(ggplot2) 
library(ggsci) 
library(hexbin)  
library(bayesplot) 
library(kableExtra)  
library(cmdstanr)  
library(webshot2) 
library(here)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(stringr)


# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

#==============================================================================
# we load our data set
#==============================================================================

data <- read.csv(here("data", "derived", "merged-data.csv"))

#==============================================================================
# we make a pairwise correlation plot of wind at time of landfall, 
# overall curvature of hurricane path, duration of hurricane and adjusted total 
# damages
#==============================================================================

# calculate duration of each hurricane in hours

data <- data %>%
  mutate(TIME = sprintf("%04s", as.integer(TIME)),
         DATE_TIME = as.POSIXct(paste(YEAR, 
                                      MONTH, 
                                      DAY, 
                                      str_sub(TIME, 1, 2), 
                                      str_sub(TIME, 3, 4)),
                                format = "%Y %m %d %H %M"))

storm_duration <- data %>%
  group_by(STORM_ID) %>%
  summarise(duration_hours = as.numeric(difftime(max(DATE_TIME), 
                                                 min(DATE_TIME), 
                                                 units = "hours")),
            .groups = "drop")

library(geosphere)

curvature <- data %>%
  arrange(STORM_ID, DATE_TIME) %>%
  group_by(STORM_ID, NAME)%>%
  mutate(
    NEXT_LAT = lead(LAT),
    NEXT_LON = lead(LON),
    BEARING = bearing(cbind(LON, LAT), cbind(NEXT_LON, NEXT_LAT)),
    # change bearings from -180 to 180 to 0 to 360 degrees
    BEARING = BEARING + 180,
    DIFF_BEARING = abs(lead(BEARING)-BEARING)) %>%
  summarise(
    total_bearing_change = sum(DIFF_BEARING, na.rm = TRUE),
    avg_landfall_wind = mean(WIND[RECORD_ID =="L"]),
    total_damage = mean(ADJ_TOTAL_DAMAGE), # damage is the same for all rows for each storm
    duration_h = as.numeric(difftime(max(DATE_TIME), min(DATE_TIME), units = "hours")),
    .groups = "drop")


    

,
    next_lon = lead(LON),
    bearing = bearing(cbind(LON, LAT), cbind(next_lon, next_lat)),
    bearing_diff = abs(bearing - lag(bearing)),
    # Normalize angle to max 180 (angles wrap around at 360)
    bearing_diff = ifelse(bearing_diff > 180, 360 - bearing_diff, bearing_diff)
  )

%>%
  summarise(
    curvature = sum(bearing_diff, na.rm = TRUE),
    avg_wind = mean(WIND[RECORD_ID =="L"]),
    total_damage = first(ADJ_TOTAL_DAMAGE),
    duration_hours = as.numeric(difftime(max(DATE_TIME), min(DATE_TIME), units = "hours")),
    .groups = "drop"
  )









