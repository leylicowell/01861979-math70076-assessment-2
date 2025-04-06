# in this script, we perform an exploratory data analysis of our cyclone
# track data from the HURDAT2 data set
# we are interested in the number of cyclones per year, how many make landfall 
# and any seasonality effects

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(ggplot2)


#==============================================================================
# we load our data set
#==============================================================================

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))


#==============================================================================
# we want to plot the number of cyclones per year from 2000-2024 
#==============================================================================

cyclones_per_year <- hurricane_data %>%
  group_by(YEAR, STORM_ID) %>%
  summarise() %>%
  group_by(YEAR) %>%
  summarise(CYCLONE_NUM = n())

# Plot
ggplot(cyclones_per_year, aes(x = YEAR, y = CYCLONE_NUM)) +
  geom_col(fill = "steelblue") +
  labs(title = "Number of Cyclones Per Year",
       x = "Year",
       y = "Number of Cyclones") +
  theme_minimal()

# save output
ggsave(here("outputs", "cyclones-per-year.pdf"), height = 4.5, width = 6, dpi=300)


#==============================================================================
# since cyclone numbers have stayed relatively constant (besides two peaks)
# over the years, we check for a monthly seasonality component 
#==============================================================================

# calculate average number of cyclones per month (across 2000-2024)
cyclones_per_month_and_year <- hurricane_data %>%
  group_by(YEAR, MONTH, STORM_ID) %>%
  summarise() %>%
  group_by(YEAR, MONTH) %>%
  summarise(CYCLONE_NUM = n())

avg_cyclones_per_month <- cyclones_per_month_and_year %>%
  group_by(MONTH) %>%
  summarise(AVG_CYCLONES = mean(CYCLONE_NUM))

# we notice some months are missing from our data set
# this is because no cyclones occurred in those months
# we now add the missing months with an avg cyclone count of zero
all_months <- data.frame(MONTH = 1:12)
avg_cyclones_per_month <- left_join(all_months, 
                                        avg_cyclones_per_month, 
                                        by = "MONTH")

# replace month numbers with month names
avg_cyclones_per_month <- avg_cyclones_per_month %>%
  mutate(AVG_CYCLONES = replace_na(AVG_CYCLONES, 0), MONTH = month.name)


# plot average number of cyclones per month
ggplot(avg_cyclones_per_month, aes(x = factor(MONTH, levels = MONTH), 
                                   y = AVG_CYCLONES)) +
  geom_col(fill = "steelblue") +
  labs(title = "Average Number of Cyclones Per Month",
       x = "Month",
       y = "Average Number of Cyclones") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45,
                                                     vjust = 1,
                                                     hjust = 1))

# we see a significant increase in the number of cyclones in July-October

# save output
ggsave(here("outputs", "avg-cyclones-per-month.pdf"), height = 4.5, width = 6, dpi=300)

#==============================================================================
# we now want to look closer into cyclone that landfall
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
ggplot(landfall_summary, aes(x = factor(CATEGORY, levels = CATEGORY), y = NUM)) +
  geom_col(fill = "steelblue") +
  labs(title = "Number of Landfalling Cyclones by Landfall Frequency, 2000-2024",
       x = "Number of Landfalls",
       y = "Number of Cyclones") +
  theme_minimal()

# save output
ggsave(here("outputs", "nbr-cyclone-landfalls.pdf"), height = 4.5, width = 6, dpi=300)

#==============================================================================
# we also want to inspect cyclone paths for cyclones that landfall
#==============================================================================

# we load our world map data
world <- map_data("world")


# we now select all the hurricanes that landfall at least once
landfall_hurricanes <- hurricane_data %>%
  filter(RECORD_ID == "L") %>%
  distinct(STORM_ID)

landfall_data <- hurricane_data %>%
  filter(STORM_ID %in% landfall_hurricanes$STORM_ID)


landfall_data <- landfall_data %>%
  mutate(TIME = sprintf("%04s", as.integer(TIME)),
    DATE_TIME = as.POSIXct(paste(YEAR, 
                                 MONTH, 
                                 DAY, 
                                 str_sub(TIME, 1, 2), 
                                 str_sub(TIME, 3, 4)),
                      format = "%Y %m %d %H %M"))

# plot map
ggplot() +
  # Base world map 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               fill = "beige", col = "#808080", linewidth = 0.2) +  
  theme(panel.background = element_rect(fill = "lightblue", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  geom_path(data = landfall_data, aes(x = LON, y = LAT, 
                                      group = STORM_ID, color = WIND), 
            size = 0.5, alpha = 0.5) +
  scale_color_gradient(low = "green", high = "red") + 
  coord_cartesian( xlim = c(-180, 160), ylim = c(-5, 70)) + 
  labs(title = "Hurricane Tracks with Wind Speed (Pacific-Centric)",
       x = "Longitude", y = "Latitude", color = "Wind Speed (knots)") +
  theme(legend.position = "right")






