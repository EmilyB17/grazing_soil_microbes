

# arguments: data = dataframe with some number of factors
# facts: character vector with the column names of the factor
# example: EmilysAwesomeFunction(veg, facts = c("Block", "Treatment", "GrazeTime"))
EmilysAwesomeFunction <- function(data, facts) { #facts is a character vector with the column names of factors
  unique.facts <- factor(facts)
  for(i in 1:length(unique.facts)) {
    levs <- levels(data[,unique.facts[i]])
    for(j in 1:length(levs)) {
      data$x <-  ifelse(data[,unique.facts[i]] %in% levs[j], 1, 0)
      names(data)[names(data) == "x"] <- paste(facts[i], levs[j], sep = "_")
    }
  }
  return(data)
}