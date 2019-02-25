

require(tidyverse)

percentDifferences <- function(df, ids) { 
  # where ids is a list of identifying columns
  
  # create new dataframe, collect data into vertical
  outDF <- df %>% 
    gather(key = "Param", value = "Value",
         names(df[which(!names(df) %in% ids)])) %>%
    # spread GrazeTime to be horizontal
    spread(key = GrazeTime, value = Value)
  # make zeros into small numbers so we can divide by "zero"
  outDF$PRE[outDF$PRE == 0] <- 0.00001
  # add columns to calculate percent difference
  outDF$diffT1 <- ((outDF$`24H` - outDF$PRE) / outDF$PRE) * 100
  outDF$diffT2 <- ((outDF$`1WK` - outDF$PRE) / outDF$PRE) * 100
  outDF$diffT3 <- ((outDF$`4WK` - outDF$PRE) / outDF$PRE) * 100
 
  # create new dataframe, remove old GrazeTime columns
  respread <- outDF %>% 
    select(-c(PRE, `24H`, `1WK`, `4WK`)) %>% 
    # make GrazeTime vertical again
    gather(key = "diffGrazeTime", value = "diff",
                 diffT1, diffT2, diffT3) %>% 
    # spread the data back to original horizontal axis
   spread(key = Param, value = diff) 
  
  # return the respread dataframe
  return(respread)
}