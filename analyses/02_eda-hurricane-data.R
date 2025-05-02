# in this script, we perform an exploratory data analysis of our cyclone
# track data from the HURDAT2 data set
# we are interested in the number of cyclone landfalls per year
# and any seasonality effects

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)
library(scales)
library(sf)
library(rnaturalearth)


#==============================================================================
# we load our data set
#==============================================================================

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))
landfall_data <- read.csv(here("data", "derived", "landfall-data.csv"))

#==============================================================================
# we want to plot the number of landfalls per year from 1900-2024 
#==============================================================================

landfalls_per_year <- landfall_data %>%
  group_by(YEAR) %>%
  summarise(CYCLONE_NUM = n(), .groups = "drop")

# Plot
p1 <- ggplot(landfalls_per_year, aes(x = YEAR, y = CYCLONE_NUM)) +
  geom_col(fill = "steelblue") +
  labs(x = "Year",
       y = "Number of Landfalls") +
  theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white", colour = NA)
  )

p1

# save output
ggsave(here("outputs", "eda-hurricane-data", "landfalls-per-year.png"), 
       plot = p1, 
       height = 4.5, 
       width = 7, 
       dpi=600)


#==============================================================================
# we check for a monthly seasonality component 
#==============================================================================

# calculate average number of cyclone landfalls per month (across 2000-2024)
landfalls_per_month_and_year <- landfall_data %>%
  group_by(YEAR, MONTH) %>%
  summarise(LANDFALL_NUM = n(), .groups = "drop")

# we notice some months are missing from our data set
# this is because no cyclones occurred in those months
# we now add the missing months with a count of zero
all_months_years <- CJ(MONTH = 1:12, YEAR = 2000:2024)

# also replace month numbers with month names
landfalls_per_month_and_year <- left_join(all_months_years, 
                                    landfalls_per_month_and_year, 
                                    by = c("MONTH", "YEAR"))

landfalls_per_month_and_year <- landfalls_per_month_and_year %>%
  mutate(LANDFALL_NUM = replace_na(LANDFALL_NUM, 0), MONTH = month.name[MONTH])

avg_landfalls_per_month <- landfalls_per_month_and_year %>%
  group_by(MONTH) %>%
  summarise(AVG_LANDFALLS = mean(LANDFALL_NUM))


# plot average number of cyclones per month
p2 <- ggplot(avg_landfalls_per_month, aes(x = factor(MONTH, levels = month.name), 
                                         y = AVG_LANDFALLS)) +
  geom_col(fill = "steelblue") +
  labs(x = "Month",
       y = "Average Number of Landfalls") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour = NA)
        )

p2
# we see a significant increase in the number of cyclones in August-October

# save output
ggsave(here("outputs", "eda-hurricane-data", "avg-landfalls-per-month.png"), 
       plot = p2,
       height = 4.5,
       width = 7, dpi=600)

#==============================================================================
# we now want to look closer into how many times cyclones landfall
#==============================================================================

# number of landfalls for each cyclone
landfall_statistics <- hurricane_data %>%
  group_by(STORM_ID) %>%
  summarise(LANDFALL_NUM = sum(RECORD_ID == "L"))

# categorise cyclones based on how many times they landfall
landfall_statistics <- landfall_statistics %>%
  mutate(CATEGORY = case_when(
    LANDFALL_NUM == 0 ~ '0',
    LANDFALL_NUM == 1 ~ '1',
    LANDFALL_NUM == 2 ~ '2',
    LANDFALL_NUM > 2  ~ ">2"))

# Count how many cyclones fall into each category
landfall_summary <- landfall_statistics %>%
  group_by(CATEGORY) %>%
  summarise(NUM = n())

# plot results
p3 <- ggplot(landfall_summary, aes(x = factor(CATEGORY, levels = CATEGORY), y = NUM)) +
  geom_col(fill = "steelblue") +
  labs(x = "Number of Landfalls per Cyclone",
       y = "Number of Cyclones") +
  theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white", colour = NA)
  )

p3

# save output
ggsave(here("outputs","eda-hurricane-data", "nbr-cyclone-landfalls.png"), 
       plot = p3, 
       height = 4.5, 
       width = 7, 
       dpi=600)

#==============================================================================
# we also want to inspect cyclone paths for cyclones that landfall 2000-2024
#==============================================================================

# we load our world map data
world <- map_data("world")


# we now select all the hurricanes that landfall at least once
landfall_hurricanes <- landfall_data %>%
  distinct(STORM_ID)

landfall_cyclones <- hurricane_data %>%
  filter(STORM_ID %in% landfall_hurricanes$STORM_ID)


landfall_cyclones <- landfall_cyclones %>%
  mutate(TIME = sprintf("%04s", as.integer(TIME)),
    DATE_TIME = as.POSIXct(paste(YEAR, 
                                 MONTH, 
                                 DAY, 
                                 str_sub(TIME, 1, 2), 
                                 str_sub(TIME, 3, 4)),
                      format = "%Y %m %d %H %M"))

# plot map
p4 <- ggplot() +
  # world map 
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group), 
               fill = "beige", 
               col = "#808080", 
               linewidth = 0.2) +  
  theme(panel.background = element_rect(fill = "#3f7ec0", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  geom_path(data = landfall_cyclones, aes(x = LON, 
                                      y = LAT, 
                                      group = STORM_ID, 
                                      color = WIND), 
            linewidth = 0.3, alpha = 0.7) +
  scale_color_gradientn(colours = c("#94d6f8", "cyan","#fbff9f","#f8e007","#fd9920","#ff1f2d", "#ee00ff"),
                        values = rescale(c(0, 33, 63, 82, 95, 112, 137)),
                        limits = c(0, 137)) + 
  coord_cartesian( xlim = c(-175, 40), ylim = c(0, 70)) + 
  labs(x = "Longitude", 
       y = "Latitude", 
       color = "Wind Speed (knots)") +
  theme(legend.position = "right")

p4

# save output
ggsave(here("outputs", "eda-hurricane-data", "landfall-map.png"), 
       plot =p4,
       height = 4, 
       width = 10, 
       dpi=600)


#==============================================================================
# we now want to examine which countries have had the most landfalls since 2000
#==============================================================================

# we count the number of landfalls per country
landfalls_per_location <- landfall_data %>%
  group_by(LOCATION) %>%
  summarise(LANDFALL_NUM = n(), .groups = "drop")

# add landfall count to countries sf object for plotting to get geometry of 
# landfall locations
countries <- ne_countries(scale = 'medium', 
                          returnclass = c("sf"), 
                          continent = c("North America", "South America"))


landfalls_per_location <- countries[, c("name_long", "geometry")] %>%
  left_join(landfalls_per_location, 
            by = c("name_long" = "LOCATION")) 

# replace NA rows with 0
landfalls_per_location <- landfalls_per_location %>%
  mutate(LANDFALL_NUM = ifelse(is.na(LANDFALL_NUM), 0, LANDFALL_NUM))

# find 10 countries with the most landfalls and save these in a dataset
top_ten_landfall_locations <- landfalls_per_location %>%
  arrange(desc(LANDFALL_NUM))

top_ten_landfall_locations <- top_ten_landfall_locations[1:10,]

# plot map of countries most at risk
p5 <- ggplot(data = landfalls_per_location) +
  geom_sf(aes(fill = LANDFALL_NUM), color = "black", linewidth = 0.125) +
  coord_sf(xlim = c(-180, -10), ylim = c(-5, 85), expand = FALSE) +
  scale_fill_gradientn(
    colors = c("#fcf9cb", "#fef67e","#f07028", "#ff6232", "#bd1d12"),  # Light red to dark red
    values = rescale(c(0, 1, 20, 50, 150))) +
  theme_minimal() +
  labs(fill = "Number of Landfalls") +
  geom_sf_text(data = top_ten_landfall_locations, 
               aes(label = name_long), 
               size = 2.5, 
               fontface = "bold", 
               color = "black",
               nudge_x = case_when(
                 top_ten_landfall_locations$name_long == "Dominican Republic" ~ 8, 
                 top_ten_landfall_locations$name_long == "Belize" ~ 3,  
                 top_ten_landfall_locations$name_long == "Cuba" ~ 2.5, 
                 top_ten_landfall_locations$name_long == "Puerto Rico" ~ 5.2,
                 top_ten_landfall_locations$name_long == "Bahamas" ~ 7,
                 top_ten_landfall_locations$name_long == "Antigua and Barbuda" ~ 9,
                 .default = 0), 
               nudge_y = case_when(
                 top_ten_landfall_locations$name_long == "Cuba" ~ 0.75,
                 top_ten_landfall_locations$name_long == "Puerto Rico" ~ 0.85,
                 top_ten_landfall_locations$name_long == "Dominican Republic" ~ 1.75,  
                 top_ten_landfall_locations$name_long == "Bahamas" ~ -1,
                 .default = 0)) +
  theme(
    axis.text = element_blank(),  
    axis.title = element_blank(),  
    panel.background = element_rect(fill = "#bde5ff"),
    panel.grid = element_blank()
  )


p5

# save output
ggsave(here("outputs", "eda-hurricane-data", "landfall-countries.png"), 
       plot =p5,
       height = 10, 
       width = 10, 
       dpi=600)


#==============================================================================
# finally we map specific landfall locations
#==============================================================================

p6 <- ggplot(data = landfalls_per_location) +
  geom_sf(fill = "beige", color = "black", linewidth = 0.125)+
  geom_point(data = landfall_data, 
             aes(x = LON, y = LAT), 
             color = "red", 
             alpha = 0.4, 
             size = 1) +
  coord_sf(xlim = c(-180, -10), ylim = c(-5, 85), expand = FALSE) +
  theme_minimal() +
  labs(fill = "Number of Landfalls") +
  geom_sf_text(data = top_ten_landfall_locations, 
               aes(label = name_long), 
               size = 2.5, 
               fontface = "bold", 
               color = "black",
               nudge_x = case_when(
                 top_ten_landfall_locations$name_long == "Dominican Republic" ~ 8, 
                 top_ten_landfall_locations$name_long == "Belize" ~ 3,  
                 top_ten_landfall_locations$name_long == "Cuba" ~ 2.5, 
                 top_ten_landfall_locations$name_long == "Puerto Rico" ~ 5.2,
                 top_ten_landfall_locations$name_long == "Bahamas" ~ 7,
                 top_ten_landfall_locations$name_long == "Antigua and Barbuda" ~ 9,
                 .default = 0), 
               nudge_y = case_when(
                 top_ten_landfall_locations$name_long == "Cuba" ~ 0.75,
                 top_ten_landfall_locations$name_long == "Puerto Rico" ~ 0.85,
                 top_ten_landfall_locations$name_long == "Dominican Republic" ~ 1.75,  
                 top_ten_landfall_locations$name_long == "Bahamas" ~ -1,
                 .default = 0)) +
  theme(
    axis.text = element_blank(),  
    axis.title = element_blank(),  
    panel.background = element_rect(fill = "#bde5ff"),
    panel.grid = element_blank()
  )

p6

# save output
ggsave(here("outputs", "eda-hurricane-data", "landfall-locations.png"), 
       plot =p6,
       height = 10, 
       width = 10, 
       dpi=600)
  
