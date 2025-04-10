# in this script, we perform an exploratory data analysis of our cyclone
# economic data from the EM_DAT data set
# note all costs are in '000 US$

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(ggplot2)
library(stringr)
library(rnaturalearth)
library(sf)


#==============================================================================
# we load our data sets
#==============================================================================

econ_data <- read.csv(here("data", "derived", "economic-data-clean.csv"))
US_data <- read.csv(here("data", "derived", "us-economic-data-clean.csv"))
merged_data <- read.csv(here("data", "derived", "merged-data.csv"))

# quick look at data sets
str(econ_data)
str(US_data)
str(merged_data)

#==============================================================================
# we want to plot the number of cyclones per country
#==============================================================================

num_cyclone_country <- econ_data %>%
  group_by(COUNTRY) %>%
  summarise(CYCLONES = n(), .groups = "drop")

# Create the plot
p1 <- ggplot(num_cyclone_country, 
             aes(x = reorder(COUNTRY, -CYCLONES), y = CYCLONES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Country", 
       y = "Number of Cyclones", 
       title = "Number of Cyclones per Country, 2000-2024") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p1

ggsave(here("outputs", "eda-economic-data", "cyclones-per-country.pdf"), 
       plot = p1,
       height = 4, 
       width = 15, 
       dpi=300)


#==============================================================================
# plot total economic damages (adjusted) per cyclone
#==============================================================================

# for plot x-axis since some hurricanes have the same names
US_data <- US_data %>%
  mutate(NAME_YEAR = paste0(NAME, " - ", YEAR))

p2 <- ggplot(US_data, aes(x = factor(NAME_YEAR, levels = NAME_YEAR), y = ADJ_TOTAL_DAMAGE)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Hurricane (Name - Year)",
       y = "Adjusted Damage (millions of US Dollars)",
       title = "Adjusted Damage per Hurricane",
       caption = "Note: Hurricane Bill (2009) is missing from the data") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2

ggsave(here("outputs", "eda-economic-data", "us-damages-per-cyclone.pdf"), 
       plot = p2,
       height = 6, 
       width = 12, 
       dpi=300)


#==============================================================================
# plot landfalls that caused damage in the US on a map
#==============================================================================

# we load our world map data
world <- map_data("world")

# we now select all hurricane landfall locations in the US from our list
landfall_data <- merged_data %>%
  filter(RECORD_ID == "L") 
us_boundary <- ne_countries(country = "United States of America", returnclass = "sf")
landfall_sf <- st_as_sf(landfall_data, coords = c("LON", "LAT"), crs = 4326)
landfall_us <- st_join(landfall_sf, us_boundary, join = st_within)

landfall_us <- as.data.frame(landfall_us %>% filter(!is.na(sovereignt)))
landfall_us <- subset(landfall_us, select = c("STORM_ID", 
                                              "NAME" ,    
                                              "DATE" ,     
                                              "TIME"  ,
                                              "RECORD_ID", 
                                              "STATUS", 
                                              "WIND", 
                                              "PRESSURE",
                                              "YEAR", 
                                              "MONTH" ,
                                              "DAY",
                                              "COUNTRY",
                                              "SUBREGION",
                                              "LOCATION",
                                              "ADJ_TOTAL_DAMAGE"))
# add back longitude and latitude coordinates
landfall_us <- left_join(landfall_us, landfall_data, by = colnames(landfall_us))



# plot map
p3 <- ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group), 
               fill = "beige", 
               col = "#808080", 
               linewidth = 0.2) +  
  theme(panel.background = element_rect(fill = "#3f7ec0", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  geom_point(data = landfall_us, aes(x = LON, 
                                     y = LAT, 
                                     color = WIND,
                                     size = ADJ_TOTAL_DAMAGE),
             alpha = 0.7) + 
  scale_color_gradient(name = "Wind Speed (knots)", low = "yellow", high = "red")+
  scale_size_continuous(name = "Adjusted Damage ($ millions)", range = c(1, 6)) +
  coord_cartesian( xlim = c(-175, -50), ylim = c(15, 55)) + 
  labs(title = "US Landfall Hurricane Locations, 2000-2024",
       x = "Longitude", 
       y = "Latitude") +
  theme(legend.position = "right")

p3

ggsave(here("outputs", "eda-economic-data", "landfall-locations.pdf"), 
       plot =p3,
       height = 4, 
       width = 10, 
       dpi=600)





