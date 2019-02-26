### FUNCTION to calculate percent differences for a time series 


require(tidyverse)

percentDifferences <- function(df, ids, timeKey, timeLevels, level1) {
  ##ARGUMENTS:
  # df = dataframe
  # ids = identifying columns (use c() not list())
  # timeKey = column with timeseries IDs
  # timeLevels = levels of the timeseries in timeKey (use c() not list())
  # level1 = dataframe w/1 column: data is the baseline timeseries
  
  # create new dataframe, collect data into vertical
  outDF <- df %>% 
    gather(key = "Param", value = "Value",
           names(df[which(!names(df) %in% ids)])) %>%
    # spread GrazeTime to be horizontal
    spread(key = timeKey, value = Value)
  # make zeros into small numbers so we can divide by "zero"
  fixZeros <- outDF[names(outDF) %in% level1] + 0.000001
  # create new vector of IDs without timeKey column
  idsNoTimeKey <- if(timeKey %in% ids) {ids[-grep(timeKey, ids)]}
  # create a new vector of TimeLevels without level1
  timeLevelsNo1 <- if(level1 %in% timeLevels) {timeLevels[-grep(level1, timeLevels)]}
  #make output DF for loop
  outDF1 <- outDF %>% select(idsNoTimeKey, Param)
  # for loop to iterate through levels and calculate % difference
  for(i in unique(which(colnames(outDF) %in% timeLevelsNo1))) {
    # create vector of the % change
    d <- ((outDF[i] - fixZeros) / fixZeros) * 100
    colnames(d) <- paste("diff", names(d), sep = "_")
    # add values to outDF1
    outDF1<- outDF1 %>% cbind(d)
  }
  # create vector of column names to search for later
  diffnames <- grep("diff", colnames(outDF1), value = TRUE)
  # create new dataframe, remove old GrazeTime columns
  respread <- outDF %>% 
    select(-timeLevels) %>% 
    # join with output from the loop
    left_join(outDF1, by = c(idsNoTimeKey, "Param")) 
  # make GrazeTime vertical again
  respread <- respread %>%  gather(key = "diffTimeSeries", value = "diff",
                                   names(respread[which(names(respread) %in% diffnames)]))
  # spread the data back to original horizontal axis
  respread1 <- respread %>% 
    spread(key = Param, value = diff) 
  # make diffTimeSeries a factor again for easier analysis
  respread1$diffTimeSeries <- factor(respread1$diffTimeSeries)
  # return the respread dataframe
  return(respread1)
}




# ---- EXAMPLE ----
# source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/2019sourcingUnivariateData_prepped.R")
## TESTING
#df <- cn18
#ids <- c("Plot", "Block", "Treatment", "GrazeTime")
#timeKey <- "GrazeTime"
#timeLevels <- c("PRE", "24H", "1WK", "4WK")
#level1 <- "PRE"

##example <- percentDifferences(df = cn18, ids = c("Plot", "Block", "Treatment", "GrazeTime"), timeKey = "GrazeTime",timeLevels = c("PRE", "24H", "1WK", "4WK"),                              level1 = "PRE")