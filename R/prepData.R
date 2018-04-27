
# Function by Shannon Albeke salbeke@uwyo.edu
# writing function to transform Emily enzyme data to be workable

# data = dataframe of enzyme data 
prepData <- function(data) {
  data$SampDate <- NULL
  data$Plot <- as.factor(data$Plot)
  data$Block <- factor(data$Block, levels = c("1", "2", "3", "4"), ordered = TRUE)
  data$Treatment <- as.factor(data$Treatment)
  data$GrazeTime <- factor(data$GrazeTime, levels = c("PRE", "24H", "1WK", "4WK"), ordered = TRUE)
  data$Sample <- as.factor(data$Sample)
  return(data)
}
