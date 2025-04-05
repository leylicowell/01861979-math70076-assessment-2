#' Reformat our hurricane data character vector into a data frame containing values of interest
#'
#' @param data A UTF-8 encoded character vector of the lines in the data file
#'
#' @returns data frame containing values of interest for our analysis
#' @export
#'

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