# In this script, we reformat the Atlantic and Pacific hurricane data from NHC, 
# clean it and save it as a .csv file

#===============================================================================
# load packages
#===============================================================================

library(dplyr)
library(magrittr)
library(here)
library(tidyr)
library(readr)
library(stringr)

#===============================================================================
# load data
#===============================================================================

#
# HURDAT2 track data
# data source: https://www.nhc.noaa.gov/data/
#

atlantic_data <- read_lines(here("data", "raw", "atlantic-data.txt"))
pacific_data <- read_lines(here("data", "raw", "pacific-data.txt"))

#===============================================================================
# reformat data sets
#===============================================================================

#-------------------------------------------------------------------------------
# we focus on the relevant columns for our analysis to simplify our data frames
# we also reformat our .txt data into a data frame object
# the function below is only used in this .R file
#-------------------------------------------------------------------------------

reformat_txt <- function(data){
  
  STORM_ID <- c()
  NAME <- c()
  DATE <- c()
  TIME <- c()
  RECORD_ID <- c()
  STATUS <- c()
  LAT <- c()
  LON <- c()
  WIND <- c()
  PRESSURE <- c()
  
  i <- 1
  while(i < length(data)){
    header_line <- str_split_1(data[i], pattern = ",")
    header_line <- gsub(" ", "", header_line) # remove whitespace
    storm_id <- header_line[1]
    name <- header_line[2]
    n_track <- header_line[3] # number of track entries
    
    for (j in 1:n_track){
      i = i + 1
      data_line <- str_split_1(data[i], pattern = ",")
      data_line <- gsub(" ", "", data_line) # remove whitespace
      STORM_ID <- c(STORM_ID, storm_id)
      NAME <- c(NAME, name)
      DATE <- c(DATE, data_line[1])
      TIME <- c(TIME,data_line[2])
      RECORD_ID <- c(RECORD_ID, data_line[3])
      STATUS <- c(STATUS, data_line[4])
      LAT <- c(LAT, data_line[5])
      LON <- c(LON, data_line[6])
      WIND <- c(WIND, data_line[7])
      PRESSURE <- c(PRESSURE, data_line[8])
    }
    i = i + 1
  }
  
  new_data <- data.frame(
    STORM_ID = STORM_ID,
    NAME = NAME,
    DATE = DATE,
    TIME = TIME,
    RECORD_ID = RECORD_ID,
    STATUS = STATUS,
    LAT = LAT,
    LON = LON,
    WIND = WIND,
    PRESSURE = PRESSURE
  )
  
  return(new_data)
}

atlantic_data <- reformat_txt(atlantic_data)
pacific_data <- reformat_txt(pacific_data)

hurricane_data <- rbind(atlantic_data, pacific_data)

#-------------------------------------------------------------------------------
# we want to analyse hurricanes from 2000 onwards
#-------------------------------------------------------------------------------

hurricane_data <- hurricane_data %>%
  filter(as.integer(str_sub(DATE, 1, 4)) >= 2000)

#-------------------------------------------------------------------------------
# reformat specific columns to get values in correct str/numeric types as they
# were loaded as character columns
#-------------------------------------------------------------------------------

hurricane_data$YEAR = as.integer(str_sub(hurricane_data$DATE, 1, 4))
hurricane_data$MONTH = as.integer(str_sub(hurricane_data$DATE, 5, 6))
hurricane_data$DAY = as.integer(str_sub(hurricane_data$DATE, 7, 8))
hurricane_data$DATE = as.Date(hurricane_data$DATE, format = '%Y%m%d')
hurricane_data$WIND = as.integer(hurricane_data$WIND)
hurricane_data$PRESSURE = as.integer(hurricane_data$PRESSURE)

# need latitude/longitude coordinates in numerical form for analysis and plotting
lat_coord <- as.numeric(sub("([0-9]+.[0-9]+).*", "\\1", hurricane_data$LAT))
hurricane_data$LAT[grep("S", hurricane_data$LAT)] = -lat_coord[grep("S", hurricane_data$LAT)]
hurricane_data$LAT[grep("N", hurricane_data$LAT)] = lat_coord[grep("N", hurricane_data$LAT)]
hurricane_data$LAT = as.numeric(hurricane_data$LAT)

lon_coord <- as.numeric(sub("([0-9]+.[0-9]+).*", "\\1", hurricane_data$LON))
hurricane_data$LON[grep("W", hurricane_data$LON)] = -lon_coord[grep("W", hurricane_data$LON)]
hurricane_data$LON[grep("E", hurricane_data$LON)] = lon_coord[grep("E", hurricane_data$LON)]
hurricane_data$LON = as.numeric(hurricane_data$LON)


#===============================================================================
# we choose to save our new data set as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#===============================================================================

write.csv(hurricane_data, 
          file = file.path(here("data", "derived"), "hurricane-data.csv"), 
          row.names = FALSE)
