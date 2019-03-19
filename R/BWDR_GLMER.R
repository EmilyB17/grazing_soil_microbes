### BWDR Data Analysis - GLMERs
require(tidyverse)
require(ggplot2)
require(emmeans)
require(MuMIn)
require(lme4)

bw <- read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/BWDR2.txt", sep = "\t", header = TRUE, stringsAsFactors =  TRUE)
bw$GrazeTime <- factor(bw$GrazeTime, ordered = TRUE, levels = c("PRE", "24H", "1WK", "4WK"))
# median is more normally distributed than mean

## Plot variability
# how much confidence do we have in the mean BWDY?

grouped <- bw %>% 
  group_by(GrazeTime, Year, Treatment) %>% 
  summarise(mean = mean(meanBWDY),
            sd = sd(sdBWDY),
            median = median(medianBWDY),
            quant2.5 = mean(quant2.5BWDY),
            quant25 = mean(quant25BWDY),
            quant50 = mean(quant50BWDY),
            quant75 = mean(quant75BWDY),
            quant97.5 = mean(quant97.5BWDY))
ggplot(data = grouped, aes(x = Treatment, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = quant2.5, ymax = quant97.5)) +
  facet_wrap(~ GrazeTime+Year)

## average by Plot
avg <- bw %>% group_by(Plot, Year, Treatment, GrazeTime, Block) %>% 
  summarise(mean = mean(meanBWDY),
            sd = sd(sdBWDY),
            median = median(medianBWDY),
            quant2.5 = mean(quant2.5BWDY),
            quant25 = mean(quant25BWDY),
            quant50 = mean(quant50BWDY),
            quant75 = mean(quant75BWDY),
            quant97.5 = mean(quant97.5BWDY))
avg$Plot <- factor(avg$Plot)
avg$Year <- factor(avg$Year)

## ---- Percent Difference ----
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")
bw18 <- avg[avg$Year %in% "2018",]
bw18 <- bw18[,c("Plot", "Treatment", "GrazeTime", "Block", "mean", "sd", "median")]
bw18.pd <- percentDifferences(df = bw18, 
                            ids = c("Plot", "Block", "Treatment", "GrazeTime"),
                            timeKey = "GrazeTime", 
                            timeLevels = c("PRE", "24H", "1WK", "4WK"),
                            level1 = "PRE")
bw18.pd$diffTimeSeries <- factor(bw18.pd$diffTimeSeries, ordered = TRUE,
                                 levels = c("diff_24H", "diff_1WK", "diff_4WK"))
bw18.pd4 <- bw18.pd[bw18.pd$Block %in% 1:3,]
bw17 <- avg[avg$Year %in% "2017",]
bw17 <- bw17[,c("Plot", "Treatment", "GrazeTime", "Block", "mean", "sd", "median")]
bw17.pd <- percentDifferences(df = bw17, 
                              ids = c("Plot", "Block", "Treatment", "GrazeTime"),
                              timeKey = "GrazeTime", 
                              timeLevels = c("PRE", "24H", "1WK", "4WK"),
                              level1 = "PRE")
bw17.pd$diffTimeSeries <- factor(bw17.pd$diffTimeSeries, ordered = TRUE,
                                 levels = c("diff_24H", "diff_1WK", "diff_4WK"))
bw17.pd4 <- bw17.pd[bw17.pd$Block %in% 1:3,]
# gather to Param column for loop
bw18.pdv <- bw18.pd4 %>% gather(key = "Param", value = "value", mean, sd, median)
bw18.pdv$log_value <- log1p(bw18.pdv$value + 101)
bw17.pdv <- bw17.pd4 %>% gather(key = "Param", value = "value", mean, sd, median)
bw17.pdv$log_value <- log1p(bw17.pdv$value + 101)

## ---2018 GLMs----

bw18.pdv$Param <- factor(bw18.pdv$Param)
params <- unique(bw18.pdv$Param)
outDF1 <- data.frame()
outDF2 <- data.frame()
for(i in 1:length(params)) {
  mod <- glmer(log_value ~ diffTimeSeries * Treatment + (1|Block),
               data = filter(bw18.pdv, Param == params[i]),
               family = gaussian(link = "log"))
  e <- as.data.frame(
    emmeans(mod, pairwise ~ Treatment | diffTimeSeries, type = "response"))
  df <- e %>% dplyr::select(c(contrasts.contrast, contrasts.diffTimeSeries,
                              contrasts.SE, contrasts.z.ratio, 
                              contrasts.p.value))
  df$Param <- factor(as.character(params[i]))
  outDF1 <- rbind(df, outDF1)
  e2 <- as.data.frame(
    emmeans(mod, pairwise ~ diffTimeSeries | Treatment, type = "response"))
  df1 <- e2 %>% dplyr::select(c(contrasts.contrast, contrasts.Treatment,
                                contrasts.SE, contrasts.z.ratio, 
                                contrasts.p.value))
  df1$Param <- factor(as.character(params[i]))
  outDF2 <- rbind(df1, outDF2)
  
}
colnames(outDF1) <- c("contrast", "contTS.Trt", "SE", "z.ratio", 
                      "p.value", "Parameter")
outDF1$sign <- ifelse(outDF1$p.value <= 0.05,
                      "significant",
                      "not_significant")
colnames(outDF2) <- c("contrast", "contTS.Trt", "SE", "z.ratio",
                      "p.value", "Parameter")
outDF2$sign <- ifelse(outDF2$p.value <= 0.05,
                      "significant",
                      "not_significant")
sig.TS1 <- outDF1[which(outDF1$sign %in% "significant"),]
sig.Trt1 <- outDF2[which(outDF2$sign %in% "significant"),]
outDF1$Year <- "2018"
outDF2$Year <- "2018"

## ----2017 GLMs----
bw17.pdv$Param <- factor(bw17.pdv$Param)
params <- unique(bw17.pdv$Param)
outDF3 <- data.frame()
outDF4 <- data.frame()
for(i in 1:length(params)) {
  mod <- glmer(log_value ~ diffTimeSeries * Treatment + (1|Block),
               data = filter(bw17.pdv, Param == params[i]),
               family = gaussian(link = "log"))
  e <- as.data.frame(
    emmeans(mod, pairwise ~ Treatment | diffTimeSeries, type = "response"))
  df <- e %>% dplyr::select(c(contrasts.contrast, contrasts.diffTimeSeries,
                              contrasts.SE, contrasts.z.ratio, 
                              contrasts.p.value))
  df$Param <- factor(as.character(params[i]))
  outDF3 <- rbind(df, outDF3)
  e2 <- as.data.frame(
    emmeans(mod, pairwise ~ diffTimeSeries | Treatment, type = "response"))
  df1 <- e2 %>% dplyr::select(c(contrasts.contrast, contrasts.Treatment,
                                contrasts.SE, contrasts.z.ratio, 
                                contrasts.p.value))
  df1$Param <- factor(as.character(params[i]))
  outDF4 <- rbind(df1, outDF4)
  
}
colnames(outDF3) <- c("contrast", "contTS.Trt", "SE", "z.ratio", 
                      "p.value", "Parameter")
outDF3$sign <- ifelse(outDF3$p.value <= 0.05,
                      "significant",
                      "not_significant")
colnames(outDF4) <- c("contrast", "contTS.Trt", "SE", "z.ratio",
                      "p.value", "Parameter")
outDF4$sign <- ifelse(outDF4$p.value <= 0.05,
                      "significant",
                      "not_significant")
sig.TS3 <- outDF3[which(outDF3$sign %in% "significant"),]
sig.Trt4 <- outDF4[which(outDF4$sign %in% "significant"),]

outDF3$Year <- "2017"
outDF4$Year <- "2017"

## ---combining both years ----
outall1 <- rbind(outDF1, outDF2, outDF3, outDF4)
# we don't care about sd so can get rid of that
gl <- outall1[!outall1$Parameter %in% "sd",]
gl <- gl[!gl$Parameter %in% "mean",]

## ---- GLM outputs: plotting significance ----
#write.table(outall1, file = "C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/BWDR_GLMER_output_NOBlock4.txt", sep = "\t", row.names = FALSE)
gl <- read.table(file = "C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/BWDR_GLMER_output_NOBlock4.txt", sep = "\t", header = TRUE)
gl <- gl[gl$Parameter %in% "median",]

gl$p.value <- round(gl$p.value, 3)
gl$CountSign <- ifelse(gl$sign %in% "significant", 1, 0)
sigs <- sum(gl$CountSign)
perc <- sum(gl$CountSign) / length(gl$CountSign) * 100 # 20% significance?

## BOXPLOTS
## PLOTTING SIGNIFICANT TREATMENTS: BOXPLOTS
sig <- gl[gl$CountSign == 1,]
sigTrt <- sig[sig$contTS.Trt %in% "diff_1WK" | 
                sig$contTS.Trt %in% "diff_24H" |
                sig$contTS.Trt %in%  "diff_4WK",]
Trtsigs <- unique(sigTrt$Parameter) 

require(ggplot2)


## 2018
unpars18 <- unique(as.character(sigTrt$Parameter[sigTrt$Parameter %in% bw18.pdv$Param & sigTrt$Year == "2018"]))
for(i in 1:length(unpars18)) {
  p <- ggplot(data = filter(bw18.pdv, Param == unpars18[i]),
              aes(x = diffTimeSeries, y = value, color = Treatment)) +
    geom_boxplot() +
    labs(x = "Time After Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("2018", unpars18[i], "BWDR", sep = " ")) +
    theme_bw()
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2018_", unpars18[i], "BWDR", "_Boxplot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}


## There is no significance for 2017

### PLOTTING SIGNIFICANT TIMES: LINEPLOTS


sigTS <- sig[sig$contTS.Trt %in% "HI" | 
               sig$contTS.Trt %in% "LO" |
               sig$contTS.Trt %in%  "NO",]
unTS18 <- unique(as.character(
  sigTS$Parameter[sigTS$Parameter %in% bw18.pdv$Param & sigTS$Year == "2018"]))
df <- bw18.pdv %>% 
  dplyr::select(-log_value) %>% 
  #spread(key = Param, value = value) %>% 
  group_by(Treatment, diffTimeSeries, Param) %>% 
  summarise(mean = mean(value),
            sd = sd(value))
dfg <- df %>%  mutate(ci.hi = round(mean + (1.96 * 
                                              (sd / sqrt(3))), 2),
                      ci.lo = round(mean - (1.96 * 
                                              (sd / sqrt(3))), 2))
for(j in 1:length(unTS18)) {
  p <- ggplot(data = filter(dfg, Param == unTS18[j]), 
              aes(x = diffTimeSeries, y = mean)) +
    geom_point(aes(color = Treatment)) +
    geom_line(aes(color = Treatment, group = Treatment)) +
    geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi, width = 0.3)) +
    facet_wrap(~ Treatment) +
    labs(x = "Time after Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("2018", unTS18[j], "BWDR", sep = " "))
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2018_", unTS18[j], "BWDR", "_LinePlot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}

## No significance for 2017