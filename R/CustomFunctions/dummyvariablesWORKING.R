

# arguments: data = dataframe with some number of factors
# facts: character vector with the column names of the factor
# example: EmilysAwesomeFunction(veg, facts = c("Block", "Treatment", "GrazeTime"))
EmilysAwesomeFunction <- function(data, facts) { 
  # how many factors are in the dataframe
  unique.facts <- unique(facts)
  # first loop to iterate through each factor
  for(i in 1:length(unique.facts)) {
    # get the levels of the specific factor
    levs <- levels(data[, unique.facts[i]])
    # second loop to iterate through each level
    for(j in 1:length(levs)) {
      # make a new column in the dataframe that stores the binary values for the level
      data$x <-  ifelse(data[,unique.facts[i]] %in% levs[j], 1, 0)
      # change the name of the column to make sense
      names(data)[names(data) == "x"] <- paste(unique.facts[i], levs[j], sep = "_")
    }
    # delete the original column so we only have the dummy variables columns
    data[, unique.facts[i]] <- NULL
  }
  # return the dataframe
  return(data)
}


# Testing!!
prepDataGLM18 <- function(data) {
  data$Plot <- NULL
  data$Block <- factor(data$Block)
  data$Treatment <- as.factor(data$Treatment)
  data$GrazeTime <- factor(data$GrazeTime)
  return(data)
}
cn <- prepDataGLM18(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018CN_data_updated10_18.txt", header = TRUE))
test <- EmilysAwesomeFunction(cn, facts = c("GrazeTime", "Block", "Treatment"))
str(test)

