# in this script, we perform an exploratory data analysis of our cyclone
# economic data from the EM_DAT data set
# note all costs are in '000 US$

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(ggplot2)
library(stringr)


#==============================================================================
# we load our data set
#==============================================================================

econ_data <- read.csv(here("data", "derived", "economic-data-clean.csv"))

# quick look at data set
str(econ_data)

#==============================================================================
# we want to plot the number of cyclones per country
#==============================================================================

num_cyclone_country <- econ_data %>%
  group_by(COUNTRY) %>%
  summarise(CYCLONES = n(), .groups = "drop")

# Create the plot
p1 <- ggplot(num_cyclone_country, aes(x = reorder(COUNTRY, -CYCLONES), y = CYCLONES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "COUNTRY", y = "Number of Cyclones", title = "Number of Cyclones per Country, 2000-2024") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p1

ggsave(here("outputs", "eda-economic-data", "cyclones-per-country.pdf"), 
       plot = p1,
       height = 4, 
       width = 15, 
       dpi=300)


#==============================================================================
# plot total economic damages (adjusted) due to cyclones per country
#==============================================================================

#------------------------------------------------------------------------------
# we notice a large number of NAs in our dataset with no indication of whether
# these represent zero values or missing data
#------------------------------------------------------------------------------



total_damage_country <- subset(econ_data, select = c())

md.pattern(econ_data)

imp <- mice(econ_data, where = is.na(econ_data$TOTAL_DAMAGE), m = 20)
total_damage_country <- econ_data %>%
  group_by(COUNTRY) %>%
  summarise(sum(!is.na(ADJ_TOTAL_DAMAGE), .groups = "drop")

# Create the plot
p1 <- ggplot(num_cyclone_country, aes(x = reorder(COUNTRY, -cyclones), y = cyclones)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Country", y = "Number of Cyclones", title = "Number of Cyclones per Country, 2000-2024") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p1

ggsave(here("outputs", "eda-economic-data", "cyclones-per-country.pdf"), 
       plot = p1,
       height = 4, 
       width = 15, 
       dpi=300)
