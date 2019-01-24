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

rpm17 <- prepData(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/rpm_Plot48.txt", sep = "\t", header = TRUE))
rpm17$biomass_kg_plot <- (82.322 * rpm17$Reading) - 341.742
rpm17[,(4:7)] <- NULL

# create microbial efficiency calculation
eff17 <- merge(mbc17, enz17, by = c("Block", "Treatment", "GrazeTime"))
eff17$efficiency <- eff17$Enzyme_nm_g_hr / eff17$MBC_mgkgdrysoil
eff17$efficiency[eff17$efficiency %in% "Inf"] <- 0
eff17$efficiency[eff17$efficiency %in% "NaN"] <- 0
eff17$Enzyme_nm_g_hr <- NULL
#eff17$efficiency[eff17$efficiency < 0] <- 0

## Forage Utlization (how much did the cattle eat of the available forage?)
# HI: PRE - 24H / PRE (24H after HI grazing)
for.hi <- data.frame(((veg$biomass_kg_plot[veg$Treatment %in% "HI" & veg$GrazeTime %in% "PRE"] -
                         veg$biomass_kg_plot[veg$Treatment %in% "HI" & veg$GrazeTime %in% "24H"]) /
                        veg$biomass_kg_plot[veg$Treatment %in% "HI" & veg$GrazeTime %in% "PRE"]) * 100)
for.hi$Treatment <- "HI"
colnames(for.hi) <- c("forage_ut", "Treatment")
# LO: PRE - 1WK / PRE (24H after LO grazing)
for.lo <- data.frame(((veg$biomass_kg_plot[veg$Treatment %in% "LO" & veg$GrazeTime %in% "PRE"] -
                         veg$biomass_kg_plot[veg$Treatment %in% "LO" & veg$GrazeTime %in% "1WK"]) /
                        veg$biomass_kg_plot[veg$Treatment %in% "LO" & veg$GrazeTime %in% "PRE"]) * 100)
for.lo$Treatment <- "LO"
colnames(for.lo) <- c("forage_ut", "Treatment")
# combine into one dataframe for comparisons
for.ut18 <- rbind(for.hi, for.lo)

## 2017
## Forage Utlization (how much did the cattle eat of the available forage?)
# HI: PRE - 24H / PRE (24H after HI grazing)
for.hi <- data.frame(((rpm17$biomass_kg_plot[rpm17$Treatment %in% "HI" & rpm17$GrazeTime %in% "PRE"] -
                         rpm17$biomass_kg_plot[rpm17$Treatment %in% "HI" & rpm17$GrazeTime %in% "24H"]) /
                        rpm17$biomass_kg_plot[rpm17$Treatment %in% "HI" & rpm17$GrazeTime %in% "PRE"]) * 100)
for.hi$Treatment <- "HI"
colnames(for.hi) <- c("forage_ut", "Treatment")
# LO: PRE - 1WK / PRE (24H after LO grazing)
for.lo <- data.frame(((rpm17$biomass_kg_plot[rpm17$Treatment %in% "LO" & rpm17$GrazeTime %in% "PRE"] -
                         rpm17$biomass_kg_plot[rpm17$Treatment %in% "LO" & rpm17$GrazeTime %in% "1WK"]) /
                        rpm17$biomass_kg_plot[rpm17$Treatment %in% "LO" & rpm17$GrazeTime %in% "PRE"]) * 100)
for.lo$Treatment <- "LO"
colnames(for.lo) <- c("forage_ut", "Treatment")
# combine into one dataframe for comparisons
for.ut17 <- rbind(for.hi, for.lo)

## 2018
## Vegetation Recovery Per Day (how much the forage grew back based on days after grazing to final sampling)
# HI: 28 recovery days
# LO: 22 recovery days
# NO: 35 recovery days (all days of the grazing trial from PRE to 4WK)

# calculate vegetation recovery
hi.rec <- data.frame((veg$biomass_kg_plot[veg$Treatment %in% "HI" & veg$GrazeTime %in% "4WK"] -
                        veg$biomass_kg_plot[veg$Treatment %in% "HI" & veg$GrazeTime %in% "24H"]) / 28)
colnames(hi.rec) <- "veg_rec"
hi.rec$Treatment <- "HI"
lo.rec <- data.frame((veg$biomass_kg_plot[veg$Treatment %in% "LO" & veg$GrazeTime %in% "4WK"] -
                        veg$biomass_kg_plot[veg$Treatment %in% "LO" & veg$GrazeTime %in% "1WK"]) / 22)
colnames(lo.rec) <- "veg_rec"
lo.rec$Treatment <- "LO"
no.rec <- data.frame((veg$biomass_kg_plot[veg$Treatment %in% "NO" & veg$GrazeTime %in% "4WK"] -
                        veg$biomass_kg_plot[veg$Treatment %in% "NO" & veg$GrazeTime %in% "PRE"]) / 35)
colnames(no.rec) <- "veg_rec"
no.rec$Treatment <- "NO"
# create data frame
veg.rec18 <- rbind(hi.rec, lo.rec, no.rec)
veg.rec18$Treatment <- as.factor(veg.rec18$Treatment)

## 2017
## Vegetation Recovery
# HI: 22 days recovery
# LO: 20 days recovery
# NO: 35 days recovery
# calculate vegetation recovery
hi.rec <- data.frame((rpm17$biomass_kg_plot[rpm17$Treatment %in% "HI" & rpm17$GrazeTime %in% "4WK"] -
                        rpm17$biomass_kg_plot[rpm17$Treatment %in% "HI" & rpm17$GrazeTime %in% "24H"]) / 22)
colnames(hi.rec) <- "veg_rec"
hi.rec$Treatment <- "HI"
lo.rec <- data.frame((rpm17$biomass_kg_plot[rpm17$Treatment %in% "LO" & rpm17$GrazeTime %in% "4WK"] -
                        rpm17$biomass_kg_plot[rpm17$Treatment %in% "LO" & rpm17$GrazeTime %in% "1WK"]) / 20)
colnames(lo.rec) <- "veg_rec"
lo.rec$Treatment <- "LO"
no.rec <- data.frame((rpm17$biomass_kg_plot[rpm17$Treatment %in% "NO" & rpm17$GrazeTime %in% "4WK"] -
                        rpm17$biomass_kg_plot[rpm17$Treatment %in% "NO" & rpm17$GrazeTime %in% "PRE"]) / 35)
colnames(no.rec) <- "veg_rec"
no.rec$Treatment <- "NO"
# create data frame
veg.rec17 <- rbind(hi.rec, lo.rec, no.rec)
veg.rec17$Treatment <- as.factor(veg.rec17$Treatment)

## Combining 2017 and 2018 data


## cn and cn17
grav18 <- cn[ , c("GrazeTime", "Block", "Treatment", "grav_mois")]
cn$grav_mois <- NULL
cn17$mineralN_mgkgdrysoil <- cn17$minN_mgkgdrysoil
cn17$minN_mgkgdrysoil <- NULL
unique(names(cn)) %in% unique(names(cn17))

# average 2017 data by plot 

cn17.avg <- data.frame(summarise(group_by(.data = cn17, Block, Treatment, GrazeTime),
                                 NO3_mgkgdrysoil = mean(NO3_mgkgdrysoil, na.rm = TRUE),
                                 NH4_mgkgdrysoil = mean(NH4_mgkgdrysoil, na.rm = TRUE),
                                 NPOC_mgkgdrysoil = mean(NPOC_mgkgdrysoil, na.rm = TRUE),
                                 DON_mgkgdrysoil = mean(DON_mgkgdrysoil, na.rm = TRUE),
                                 mineralN_mgkgdrysoil = mean(mineralN_mgkgdrysoil, na.rm = TRUE)))
cn.all <- rbind(cn17.avg, cn)

## Gravimetric Moisture
grav17.avg <- data.frame(summarise(group_by(.data = grav17, Block, Treatment, GrazeTime),
                                   grav_mois = mean(Gravmois, na.rm = TRUE)))
grav.all <- rbind(grav18, grav17.avg)

## Enzymes
enz17.avg <- data.frame(summarise(group_by(.data = enz17, Block, Treatment,
                                           GrazeTime, Substrate),
                                  Enzyme_nm_g_hr = mean(Enzyme_nm_g_hr, na.rm = TRUE)))
enz.all <- rbind(enz17.avg, enz)

## Veg Biomass
veg$reading_rpm <- NULL
veg17.avg <- data.frame(summarise(group_by(.data = rpm17, Block, Treatment, GrazeTime),
                                  biomass_kg_plot = mean(biomass_kg_plot, na.rm = TRUE)))
veg.all <- rbind(veg, veg17.avg)

## Forage Utilization
for.ut.all <- rbind(for.ut17, for.ut18)

## Vegetation Recovery
veg.rec.all <- rbind(veg.rec17, veg.rec18)