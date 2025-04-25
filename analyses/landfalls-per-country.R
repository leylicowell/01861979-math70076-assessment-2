# in this script, we 

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(stringr)
library(data.table)
library(rnaturalearth)
library(sf)

#==============================================================================
# we load our data set
#==============================================================================

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))

#==============================================================================
# we add countries to our landfall data
#==============================================================================

landfalls <- hurricane_data %>%
  filter(RECORD_ID == "L")

landfalls_sf <- st_as_sf(landfalls, coords = c("LON", "LAT"), crs = 4326)

countries <- ne_countries(scale = 'medium', returnclass = c("sf")) 

landfall_countries <- st_join(landfalls_sf, countries)

landfalls$COUNTRY <- landfall_countries$name_long

#-------------------------------------------------------------------------------
# landfall is defined as the center of the system crossing a coastline 
# according to the NOAA, the hurricane's center usually spans 32-64 km
# for this reason, some of our landfall locations are assigned NA values as they 
# don't lie within country boundaries
# for these data points, we assign the nearest country to the landfall location
#-------------------------------------------------------------------------------

missing_idx <- which(is.na(landfalls$COUNTRY))
nearest_countries_idx <- st_nearest_feature(landfalls_sf[missing_idx, ], countries)
landfalls$COUNTRY[missing_idx] <- countries$name_long[nearest_countries_idx]

#==============================================================================
# save data set
#==============================================================================

write.csv(landfalls, 
          file = file.path(here("data", "derived"), "landfall-data.csv"), 
          row.names = FALSE)

