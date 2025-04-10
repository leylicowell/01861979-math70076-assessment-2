# In this script, we load and clean the economic data from EM-DAT, 2000-2024,
# we merge it with our data of hurricanes with billion dollar damages as this 
# contains values missing from our initial EM-DAT data set, 
# and save it as a .csv file
# we also save a separate data set containing US damages for our statistical models

#------------------------------------------------------------------------------
# load packages and functions
#------------------------------------------------------------------------------

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(readr)
library(stringr)
library(readxl)

#------------------------------------------------------------------------------
# load data
#------------------------------------------------------------------------------

#
# EM-DAT tropical cyclones in the Americas 
# data source: https://www.emdat.be
#

econ_data <- read_xlsx(here("data", "raw", "economic-data.xlsx"))
US_econ_data <- read_csv(here("data", "raw", "us-economic-data.csv"), skip = 2)

#------------------------------------------------------------------------------
# we focus on the relevant columns for our analysis to simplify our data frames
# note this data set only contains data from 2000 onwards
#------------------------------------------------------------------------------

# standardise column names for readability
colnames(econ_data) <- toupper(colnames(econ_data))
colnames(US_econ_data) <- toupper(colnames(US_econ_data))

# select columns of interest only
econ_data <- econ_data[c("DISNO.", 
                         "EVENT NAME",
                         "COUNTRY", 
                         "SUBREGION", 
                         "LOCATION", 
                         "START YEAR",
                         "TOTAL DAMAGE, ADJUSTED ('000 US$)")]

US_econ_data <- US_econ_data[c("NAME", "BEGIN DATE", "CPI-ADJUSTED COST")]


# reformat column names
colnames(econ_data) = c("ID", 
                        "NAME",
                        "COUNTRY", 
                        "SUBREGION", 
                        "LOCATION", 
                        "YEAR",
                        "ADJ_TOTAL_DAMAGE")

colnames(US_econ_data) = c("NAME", "YEAR", "ADJ_TOTAL_DAMAGE", "TOTAL_DAMAGE")

# econ_data values are in thousands of US$ so we convert the appropriate columns to millions 
econ_data <- econ_data %>%
  mutate(ADJ_TOTAL_DAMAGE = ADJ_TOTAL_DAMAGE/1e3)

# START_YEAR in US_econ_data is currently a date so we want to extract the year
US_econ_data <- US_econ_data %>%
  mutate(START_YEAR = as.integer(str_sub(as.character(YEAR),1,4)))


#------------------------------------------------------------------------------
# US specific
#------------------------------------------------------------------------------

US_data <- econ_data %>%
  filter(COUNTRY == "United States of America")

# we want to retrieve the US hurricane names for which we have NA values 
# in total adjusted damage from US_econ_data but to do this we first need to
# standardise names in both datasets

US_econ_data$NAME <- gsub("\\(.*\\)", "", US_econ_data$NAME)
US_econ_data$NAME <- gsub("Hurricane","", US_econ_data$NAME)
US_econ_data$NAME <- gsub("Tropical Storm","", US_econ_data$NAME)
US_econ_data$NAME <- gsub("Typhoon","", US_econ_data$NAME)
US_econ_data$NAME <- gsub(" ", "", US_econ_data$NAME)

US_data$NAME <- gsub(".*'(.*?'.*?').*", "\\1",US_data$NAME)
US_data$NAME <- gsub(".*'(.*?)'.*", "\\1",US_data$NAME)
US_data$NAME <- gsub(".*\"(.*?)\".*", "\\1", US_data$NAME)
US_data$NAME <- gsub("Hurricane", "\\1", US_data$NAME)
US_data$NAME <- gsub("Tropical storm", "\\1", US_data$NAME)
US_data$NAME <- gsub(" ", "", US_data$NAME)


# we retrieve the names of the hurricanes with missing ADJ_TOTAL_DAMAGE values
missing_hurricanes <- US_data[which(is.na(US_data$ADJ_TOTAL_DAMAGE)),]
print(missing_hurricanes$NAME)

damage_values <- US_econ_data[US_econ_data$NAME %in% missing_hurricanes$NAME, ]$ADJ_TOTAL_DAMAGE

US_data$ADJ_TOTAL_DAMAGE[which(is.na(US_data$ADJ_TOTAL_DAMAGE) & US_data$NAME != "Bill")] = damage_values

# no info found online regarding hurricane Bill damages and since data is MNAR 
# standard imputation methods lead to high bias
# multiple imputation also leads to significantly incorrect results (damages in billions)
# which we know to be incorrect due to its absence from our US_econ_data 
# (which is a list of all hurricanes causing >1 billion $ in damages)
# as such, since this is only one data point, we choose to ignore it in our data for now,
# and to try and estimate it using a bayesian regression model we will fit in our analysis

#------------------------------------------------------------------------------
# we choose to save our new data set as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#------------------------------------------------------------------------------

write.csv(econ_data, 
          file = file.path(here("data", "derived"), "economic-data-clean.csv"), 
          row.names = FALSE)


write.csv(US_data, 
          file = file.path(here("data", "derived"), "us-economic-data-clean.csv"), 
          row.names = FALSE)
