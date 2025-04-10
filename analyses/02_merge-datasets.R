# in this script, we merge our US economic data with our HURDAT2 data set

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(stringr)
library(data.table)
library(mice)

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

# merge datasets
merged_data <- merge(hurricane_data, us_econ_data, by = c("NAME", "YEAR"))

#------------------------------------------------------------------------------
# we choose to save our new data set as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#------------------------------------------------------------------------------

write.csv(merged_data, 
          file = file.path(here("data", "derived"), "merged-data.csv"), 
          row.names = FALSE)