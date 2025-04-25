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


#==============================================================================
# we load our data set
#==============================================================================

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))

#==============================================================================
# we want to plot the number of landfalls per year from 1900-2024 
#==============================================================================

landfalls_per_year <- hurricane_data %>%
  filter(RECORD_ID == "L") %>%
  group_by(YEAR) %>%
  summarise(CYCLONE_NUM = n(), .groups = "drop")

# Plot
p1 <- ggplot(landfalls_per_year, aes(x = YEAR, y = CYCLONE_NUM)) +
  geom_col(fill = "steelblue") +
  labs(title = "Number of Cyclone Landfalls Per Year",
       x = "Year",
       y = "Number of Landfalls") +
  theme_minimal()

p1

# save output
ggsave(here("outputs", "eda-hurricane-data", "landfalls-per-year.pdf"), 
       plot = p1, 
       height = 4.5, 
       width = 7, 
       dpi=300)


#==============================================================================
# we check for a monthly seasonality component 
#==============================================================================

# calculate average number of cyclone landfalls per month (across 2000-2024)
landfalls_per_month_and_year <- hurricane_data %>%
  filter(RECORD_ID == "L") %>%
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
  labs(title = "Average Number of Landfalls Per Month, averaging across 2000-2024",
       x = "Month",
       y = "Average Number of Landfalls") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))

p2
# we see a significant increase in the number of cyclones in August-October

# save output
ggsave(here("outputs", "eda-hurricane-data", "avg-landfalls-per-month.pdf"), 
       plot = p2,
       height = 4.5,
       width = 7, dpi=300)

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
  labs(title = "Number of Landfalling Cyclones by Landfall Frequency, 2000-2024",
       x = "Number of Landfalls per Cyclone",
       y = "Number of Cyclones") +
  theme_minimal()

p3
# save output
ggsave(here("outputs","eda-hurricane-data", "nbr-cyclone-landfalls.pdf"), 
       plot = p3, 
       height = 4.5, 
       width = 7, 
       dpi=300)

#==============================================================================
# we also want to inspect cyclone paths for cyclones that landfall 2000-2024
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
p4 <- ggplot() +
  # world map 
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group), 
               fill = "beige", 
               col = "#808080", 
               linewidth = 0.2) +  
  theme(panel.background = element_rect(fill = "#3f7ec0", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  geom_path(data = landfall_data, aes(x = LON, 
                                      y = LAT, 
                                      group = STORM_ID, 
                                      color = WIND), 
            linewidth = 0.3, alpha = 0.7) +
  scale_color_gradientn(colours = c("#94d6f8", "cyan","#fbff9f","#f8e007","#fd9920","#ff1f2d", "#ee00ff"),
                        values = rescale(c(0, 33, 63, 82, 95, 112, 137)),
                        limits = c(0, 137)) + 
  coord_cartesian( xlim = c(-175, 40), ylim = c(0, 70)) + 
  labs(title = "Landfall Hurricane Tracks with Wind Speed, 2000-2024",
       x = "Longitude", 
       y = "Latitude", 
       color = "Wind Speed (knots)") +
  theme(legend.position = "right")

p4

ggsave(here("outputs", "eda-hurricane-data", "landfall-map.pdf"), 
       plot =p4,
       height = 4, 
       width = 10, 
       dpi=600)



