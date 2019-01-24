## All packages required & sourced data
if(!require(prettydoc)) install.packages("prettydoc")
library(prettydoc)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
if(!require(lme4)) install.packages("lme4")
library(lme4)
if(!require(MuMIn)) install.packages("MuMIn")
library(MuMIn)

# prepData
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/prepData.R")
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/prepData2018.R")



## Read in data
cn18 <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018CN_data_updated10_18.txt", header = TRUE))

enz18 <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018enzymes_vertical.txt", header = TRUE))

veg18 <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018rpm.txt", header = TRUE))

cn17 <- prepData(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2017CN_data_updated10_18.txt", sep = "\t", header = TRUE))

mbc17 <- cn17[complete.cases(cn17),] # removes rows with NA values for MBC & MBN

cn17 <- cn17[, -(9:10)] # removes MBC and MBN columns so NAs don't mess up GLMs

grav17 <- prepData(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2017Gravimetric_moisture.txt", sep = "\t", header = TRUE))

enz17 <- prepData(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2017enzymes_vertical.txt", sep = "\t", header = TRUE))

ph17 <- prepData(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2017pH.txt", sep = "\t", header = TRUE))

rpm17 <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/rpm_Plot48.txt", sep = "\t", header = TRUE))
rpm17$biomass_kg_plot <- (82.322 * rpm17$Reading) - 341.742
rpm17[,(4:7)] <- NULL

# create microbial efficiency calculation
eff17 <- merge(mbc17, enz17, by = c("Block", "Treatment", "GrazeTime"))
eff17$efficiency <- eff17$Enzyme_nm_g_hr / eff17$MBC_mgkgdrysoil
eff17$efficiency[eff17$efficiency %in% "Inf"] <- 0
eff17$efficiency[eff17$efficiency %in% "NaN"] <- 0
eff17$Enzyme_nm_g_hr <- NULL
#eff17$efficiency[eff17$efficiency < 0] <- 0

