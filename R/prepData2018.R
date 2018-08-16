# Function by Emily Bean ebean@uwyo.edu and Shannon Albeke salbeke@uwyo.edu
# Writing function to transform Emily field data from 2018 to be workable

# data = dataframe of field data that must contain the following variables:
# Plot, Block, Treatment, GrazeTime

prepData2018 <- function(data) {
  data$Plot <- as.factor(data$Plot)
  data$Block <- factor(data$Block, levels = c("1", "2", "3", "4"), ordered = TRUE)
  data$Treatment <- as.factor(data$Treatment)
  data$GrazeTime <- factor(data$GrazeTime, levels = c("PRE", "24H", "1WK", "4WK"), ordered = TRUE)
  return(data)
}

