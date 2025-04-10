# in this script, we merge our US economic data with our HURDAT2 data set

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
us_econ_data <- read.csv(here("data", "derived", "us-economic-data-clean.csv"))

#==============================================================================
# merge data sets
#==============================================================================

# reformat us_econ_data NAME column so names are in capitals (like hurricane_data)
us_econ_data$NAME = toupper(us_econ_data$NAME)

# we now select all cyclones that landfall in the US
landfall_data <- hurricane_data %>%
  filter(RECORD_ID == "L") 
us_boundary <- ne_countries(country = "United States of America", returnclass = "sf")
landfall_sf <- st_as_sf(landfall_data, coords = c("LON", "LAT"), crs = 4326)
landfall_us <- st_join(landfall_sf, us_boundary, join = st_within)

landfall_us <- as.data.frame(landfall_us %>% filter(!is.na(sovereignt)))

# retrieve storm IDs 
storm_ids <- unique(landfall_us$STORM_ID)

# select the track data for these cyclones only
us_track_data <- hurricane_data %>%
  filter(STORM_ID %in% storm_ids)

# merge datasets
merged_data <- merge(us_track_data, us_econ_data, by = c("NAME", "YEAR"))

#------------------------------------------------------------------------------
# we choose to save our new data set as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#------------------------------------------------------------------------------

write.csv(merged_data, 
          file = file.path(here("data", "derived"), "merged-data.csv"), 
          row.names = FALSE)