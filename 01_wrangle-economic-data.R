# In this script, we load and clean the economic data from EM-DAT
# and save it as a .csv file

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

#------------------------------------------------------------------------------
# we focus on the relevant columns for our analysis to simplify our data frames
# note this data set only contains data from 2000 onwards
#------------------------------------------------------------------------------

# standardise column names for readability
colnames(econ_data) <- toupper(colnames(econ_data))
econ_data <- econ_data[c("DISNO.", "EVENT NAME","COUNTRY", "SUBREGION", 
                         "LOCATION", "ASSOCIATED TYPES", "START YEAR",
                         "START MONTH","START DAY", "END YEAR", 
                         "RECONSTRUCTION COSTS ('000 US$)", 
                         "RECONSTRUCTION COSTS, ADJUSTED ('000 US$)", 
                         "INSURED DAMAGE ('000 US$)" ,
                         "INSURED DAMAGE, ADJUSTED ('000 US$)",
                         "TOTAL DAMAGE ('000 US$)",
                         "TOTAL DAMAGE, ADJUSTED ('000 US$)","CPI")]


#------------------------------------------------------------------------------
# we choose to save our new data set as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#------------------------------------------------------------------------------

write.csv(econ_data, 
          file = file.path(here("data", "derived"), "economic-data-clean.csv"), 
          row.names = FALSE)
