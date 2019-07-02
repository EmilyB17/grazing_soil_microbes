## Indentity-linked Gaussian on log-transformed Percent Differences data

source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/2019sourcingUnivariateData_prepped.R")
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")
require(tidyverse)
require(emmeans)

all18 <- cn18 %>% 
  merge(veg18, by = c("Plot", "GrazeTime", "Block", "Treatment"))

all18$reading_rpm <- NULL
all18[,5:11] <- all18[,5:11] + 0.00001
all18.pd <- percentDifferences(df = all18,
                               ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                               timeKey = "GrazeTime",
                               timeLevels = c("PRE", "24H", "1WK", "4WK"),
                               level1 = "PRE")
all18.pd <- all18.pd[all18.pd$Block %in% 1:3,]
all18.pdv <- all18.pd %>% 
  gather(key = "Param", value = "value",
         biomass_kg_plot, DON_mgkgdrysoil, grav_mois,
         mineralN_mgkgdrysoil, NH4_mgkgdrysoil, NO3_mgkgdrysoil, NPOC_mgkgdrysoil)
all18.pdv$log_value <- log1p(all18.pdv$value + 101)

## for results section: how much did vegetation decrease by the end of grazing for each treatment?
grz <- all18.pdv %>% 
  spread(key = Param, value = value) %>% 
  group_by(diffTimeSeries, Treatment) %>% 
  summarise(meanBio = mean(biomass_kg_plot))

# ----2018 GLMs----
## GLM LOOP TO STORE MODEL PARAMETERS - 2017 Enzymes


all18.pdv$Param <- factor(all18.pdv$Param)
params <- unique(all18.pdv$Param)
outDF1 <- data.frame()
outDF2 <- data.frame()
modelFit <- data.frame()
for(i in 1:length(params)) {
  mod <- glm(log_value ~ diffTimeSeries * Treatment,
               data = filter(all18.pdv, Param == params[i]),
               family = gaussian(link = "identity"))
  mf <- data.frame(Param = as.character(params[i]),
                   Year = "2018",
                   deviance = mod$deviance,
                   null.deviance = mod$null.deviance,
                   diff = mod$null.deviance - mod$deviance,
                   df.null = mod$df.null,
                   df.dev = mod$df.residual)
  modelFit <- rbind(mf, modelFit)
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

# ---- 2017 Data Wrangling ----
myMean<- function(df, keep, kill) {
  vert <- df %>% 
    dplyr::select(-kill) 
  vert <- vert %>% # remove the columns that you don't want 
    gather(key = "Param", value = "Value", names(vert)[which(!names(vert) %in% keep)]) %>% # gather by columns you do want
    group_by(Plot, Block, GrazeTime, Treatment, Param) %>%  # group by columns you want to keep
    summarize(MeanValue = mean(Value, na.rm = TRUE)) %>% # summarize by average
    spread(key = Param, value = MeanValue) # make horizontal again
  return(vert)
}
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/prepData2018.R")
rpm17 <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/rpm_Plot48.txt", sep = "\t", header = TRUE))
rpm17$biomass_kg_plot <- (82.322 * rpm17$Reading) - 341.742
rpm17<- rpm17[,c("Plot", "Block", "Treatment", "GrazeTime", "biomass_kg_plot")]
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/prepData.R")
cn17 <- prepData(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2017CN_data_updated10_18.txt", sep = "\t", header = TRUE))


all17 <- merge(myMean(df = cn17,
                keep = c("Plot", "Block", "GrazeTime", "Treatment"),
                kill = c("Sample")),
               myMean(df = grav17,
                      keep = c("Plot", "Block", "GrazeTime", "Treatment"),
                      kill = c("Sample")),
               by = c("Plot", "Block", "GrazeTime", "Treatment")) %>% 
  merge(rpm17, by = c("Plot", "Block", "GrazeTime", "Treatment"))



# calculate % change
all17[,5:11] <- all17[,5:11] + 0.00001
all17.pd <- percentDifferences(df = all17,
                               ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                               timeKey = "GrazeTime",
                               timeLevels = c("PRE", "24H", "1WK", "4WK"),
                               level1 = "PRE")
all17.pd <- all17.pd[all17.pd$Block %in% 1:3,]
all17.pdv <- all17.pd %>% 
  gather(key = "Param", value = "value",
         biomass_kg_plot, DON_mgkgdrysoil, Gravmois,
         MBC_mgkgdrysoil, MBN_mgkgdrysoil, minN_mgkgdrysoil,
         NO3_mgkgdrysoil, NH4_mgkgdrysoil, NPOC_mgkgdrysoil)


## for results section: how much did vegetation decrease by the end of grazing for each treatment?
grz <- all17.pdv %>% 
  spread(key = Param, value = value) %>% 
  group_by(diffTimeSeries, Treatment) %>% 
  summarise(meanBio = mean(biomass_kg_plot))


all17.pdv$log_value <- log1p(all17.pdv$value + 101)
# ---- 2017 GLMs ----

all17.pdv$Param <- factor(all17.pdv$Param)
params <- unique(all17.pdv$Param)
outDF3 <- data.frame()
outDF4 <- data.frame()
for(i in 1:length(params)) {
  mod <- glm(log_value ~ diffTimeSeries * Treatment,
               data = filter(all17.pdv, Param == params[i]),
               family = gaussian(link = "identity"))
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

# ---- Combining both to write to table ----

outall1 <- rbind(outDF1, outDF2, outDF3, outDF4)



# adding Enzyme GLMERs -- source 2019Enzyme_DataAnalysis

colnames(eoutDF1) <- c("contrast", "contTS.Trt", "SE", "z.ratio",
                       "p.value", "Parameter", "sign", "Year")
colnames(eoutDF2) <- c("contrast", "contTS.Trt", "SE", "z.ratio",
                       "p.value", "Parameter", "sign", "Year")
colnames(eoutDF3) <- c("contrast", "contTS.Trt", "SE", "z.ratio",
                       "p.value", "Parameter", "sign", "Year")

colnames(eoutDF4) <- c("contrast", "contTS.Trt", "SE", "z.ratio",
                       "p.value", "Parameter", "sign", "Year")
enzglm <- rbind(eoutDF1, eoutDF2, eoutDF3, eoutDF4)

#enzglm1 <- enzglm %>% 
  #mutate(Parameter = Enzyme) %>% 
  #dplyr::select(-Enzyme)

#nm <- which(names(enzglm1) %in% names(outall))
outall.enz <- rbind(enzglm, outall1)

## WRITE TABLE
#write.table(outall.enz, file = "C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/v3_IDlink_GLMoutput2017_2018.txt", sep = "\t", row.names = FALSE)


# ----- GLMER outputs ----

gl <- read.table("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/v2GLMoutput2017_2018.txt", sep = "\t", header = TRUE)
gl$p.value <- round(gl$p.value, 3)
gl$CountSign <- ifelse(gl$sign %in% "significant", 1, 0)
sigs <- sum(gl$CountSign)
perc <- sum(gl$CountSign) / length(gl$CountSign) * 100

gleea <- read.table("//petalibrary.arcc.uwyo.edu/homes/lvandiep/SoilEcologyLab/Students/Bean/THESIS WRITING/Ch3_Soils/updatedEEAGLMS_noblock4.txt", sep = "\t", header = TRUE)
sig <- gleea[gleea$sign %in% "significant",]

## PLOTTING SIGNIFICANT TREATMENTS: BOXPLOTS
sig <- gl[gl$CountSign == 1,]
sigTrt <- sig[sig$contTS.Trt %in% "diff_1WK" | 
                sig$contTS.Trt %in% "diff_24H" |
                sig$contTS.Trt %in%  "diff_4WK",]
Trtsigs <- unique(sigTrt$Parameter) # 26 parameters have significance

require(ggplot2)
# data frames for non-enzymes: all17.pdv and all18.pdv
# already have plots for significant enzymes

## 2018
all18.pdv$diffTimeSeries <- factor(all18.pdv$diffTimeSeries, ordered = T,
                                   levels = c("diff_24H", "diff_1WK", "diff_4WK"))

unpars18 <- unique(as.character(sigTrt$Parameter[sigTrt$Parameter %in% all18.pdv$Param & sigTrt$Year == "2018"]))
for(i in 1:length(unpars18)) {
  p <- ggplot(data = filter(all18.pdv, Param == unpars18[i]),
              aes(x = diffTimeSeries, y = value, color = Treatment)) +
    geom_boxplot() +
    labs(x = "Time After Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("GLM ID link: 2018", unpars18[i], sep = " ")) +
    theme_bw()
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2018_GLMIDlink_", unpars18[i], "_Boxplot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}

## 2017
all17.pdv$diffTimeSeries <- factor(all17.pdv$diffTimeSeries, ordered = T,
                                   levels = c("diff_24H", "diff_1WK", "diff_4WK"))

unpars17 <- unique(as.character(sigTrt$Parameter[sigTrt$Parameter %in% all17.pdv$Param & sigTrt$Year == "2017"]))
for(i in 1:length(unpars17)) {
  p <- ggplot(data = filter(all17.pdv, Param == unpars17[i]),
              aes(x = diffTimeSeries, y = value, color = Treatment)) +
    geom_boxplot() +
    labs(x = "Time After Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("GLM ID link: 2017", unpars17[i], sep = " ")) +
    theme_bw()
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2017_GLMIDlink_", unpars17[i], "_Boxplot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}

### PLOTTING SIGNIFICANT TIMES: LINEPLOTS
# already have line plots for enzymes

sigTS <- sig[sig$contTS.Trt %in% "HI" | 
                sig$contTS.Trt %in% "LO" |
                sig$contTS.Trt %in%  "NO",]
unTS18 <- unique(as.character(
  sigTS$Substrate[sigTS$Substrate %in% enz2018.pd.cy$Substrate & sigTS$Year == "2018"]))
df <- enz2018.pd.cy %>% 
  dplyr::select(-log_Enz) %>% 
  #spread(key = Param, value = value) %>% 
  group_by(Treatment, diffTimeSeries, Substrate) %>% 
  summarise(mean = mean(Enzyme_nm_g_hr),
            sd = sd(Enzyme_nm_g_hr))
dfg <- df %>%  mutate(ci.hi = round(mean + (1.96 * 
                                              (sd / sqrt(3))), 2),
                      ci.lo = round(mean - (1.96 * 
                                              (sd / sqrt(3))), 2))
dfg$diffTimeSeries <- factor(dfg$diffTimeSeries, ordered = TRUE, levels = c("diff_24H", "diff_1WK", "diff_4WK"))
for(j in 1:length(unTS18)) {
  p <- ggplot(data = filter(dfg, Substrate == unTS18[j]), 
              aes(x = diffTimeSeries, y = mean)) +
    geom_point(aes(color = Treatment)) +
    geom_line(aes(color = Treatment, group = Treatment)) +
    geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi, width = 0.3)) +
    facet_wrap(~ Treatment) +
    labs(x = "Time after Grazing", y = "Percent Change from PRE") +
    ggtitle(paste(" 2018", unTS18[j], sep = " "))
  ggsave(filename = paste("//petalibrary.arcc.uwyo.edu/homes/lvandiep/SoilEcologyLab/Students/Bean/THESIS WRITING/Ch3_Soils/Figures_rough/", "EEAGLMS_", unTS18[j], "_LinePlot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}

# 2017

unTS17 <- unique(as.character(
  sigTS$Parameter[sigTS$Parameter %in% all17.pdv$Param & sigTS$Year == "2017"]))
df <- all17.pdv %>% 
  dplyr::select(-log_value) %>% 
  #spread(key = Param, value = value) %>% 
  group_by(Treatment, diffTimeSeries, Param) %>% 
  summarise(mean = mean(value),
            sd = sd(value))
dfg <- df %>%  mutate(ci.hi = round(mean + (1.96 * 
                                              (sd / sqrt(3))), 2),
                      ci.lo = round(mean - (1.96 * 
                                              (sd / sqrt(3))), 2))
dfg$diffTimeSeries <- factor(dfg$diffTimeSeries, ordered = TRUE, levels = c("diff_24H", "diff_1WK", "diff_4WK"))
for(j in 1:length(unTS17)) {
  p <- ggplot(data = filter(dfg, Param == unTS17[j]), 
              aes(x = diffTimeSeries, y = mean)) +
    geom_point(aes(color = Treatment)) +
    geom_line(aes(color = Treatment, group = Treatment)) +
    geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi, width = 0.3)) +
    facet_wrap(~ Treatment) +
    labs(x = "Time after Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("GLM ID link: 2017", unTS17[j], sep = " "))
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2017_GLMIDlink_", unTS17[j], "_LinePlotCORRECTED.jpeg", sep = ""), plot = p,
         device = "jpeg")
}


## PLOTTING ENZYME SIGNIFICANCE
#2018
enz18sigTrt <- unique(as.character(
  sigTrt$Parameter[sigTrt$Parameter %in% enz2018.pd.cy$Substrate & 
                     sigTrt$Year == "2018"]))
for(i in 1:length(enz18sigTrt)) {
  p <- ggplot(data = filter(enz2018.pd.cy, Substrate == enz18sigTrt[i]),
              aes(x = diffTimeSeries, y = Enzyme_nm_g_hr, color = Treatment))+
    geom_boxplot() +
    labs(x = "Time after Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("GLM ID link: 2018", enz18sigTrt[i], sep = " "))
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2018_GLMIDlink_", enz18sigTrt[i], "_Boxplot.jpeg", sep = ""), plot = p,
         device = "jpeg")
    
}

# 2017
enz17sigTrt <- unique(as.character(
  sigTrt$Parameter[sigTrt$Parameter %in% enz2017.pd4$Substrate & 
                     sigTrt$Year == "2017"])) 
# THERE IS NO SIGNIFICANCE


## LINEPLOTS
# 2018
sigTS <- sig[sig$contTS.Trt %in% "HI" | 
               sig$contTS.Trt %in% "LO" |
               sig$contTS.Trt %in%  "NO",]
enzTS18 <- unique(as.character(
  sigTS$Parameter[sigTS$Parameter %in% enz2018.pd.cy$Substrate & sigTS$Year == "2018"]))
df <- enz2018.pd.cy %>% 
  dplyr::select(-log_Enz) %>% 
  #spread(key = Param, value = value) %>% 
  group_by(Treatment, diffTimeSeries, Substrate) %>% 
  summarise(mean = mean(Enzyme_nm_g_hr),
            sd = sd(Enzyme_nm_g_hr))
dfg <- df %>%  mutate(ci.hi = round(mean + (1.96 * 
                                              (sd / sqrt(3))), 2),
                      ci.lo = round(mean - (1.96 * 
                                              (sd / sqrt(3))), 2))
for(j in 1:length(enzTS18)) {
  p <- ggplot(data = filter(dfg, Substrate == enzTS18[j]), 
              aes(x = diffTimeSeries, y = mean)) +
    geom_point(aes(color = Treatment)) +
    geom_line(aes(color = Treatment, group = Treatment)) +
    geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi, width = 0.3)) +
    facet_wrap(~ Treatment) +
    labs(x = "Time after Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("GLM ID link: 2018", enzTS18[j], sep = " "))
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2018_GLMIDlink_", enzTS18[j], "_LinePlot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}

# 2017
enzTS17 <- unique(as.character(
  sigTS$Parameter[sigTS$Parameter %in% enz2017.pd4$Substrate & sigTS$Year == "2017"]))
df <- enz2017.pd4 %>% 
  dplyr::select(-log_Enz) %>% 
  #spread(key = Param, value = value) %>% 
  group_by(Treatment, diffTimeSeries, Substrate) %>% 
  summarise(mean = mean(Enzyme_nm_g_hr),
            sd = sd(Enzyme_nm_g_hr))
dfg <- df %>%  mutate(ci.hi = round(mean + (1.96 * 
                                              (sd / sqrt(3))), 2),
                      ci.lo = round(mean - (1.96 * 
                                              (sd / sqrt(3))), 2))
for(j in 1:length(enzTS17)) {
  p <- ggplot(data = filter(dfg, Substrate == enzTS17[j]), 
              aes(x = diffTimeSeries, y = mean)) +
    geom_point(aes(color = Treatment)) +
    geom_line(aes(color = Treatment, group = Treatment)) +
    geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi, width = 0.3)) +
    facet_wrap(~ Treatment) +
    labs(x = "Time after Grazing", y = "Percent Change from PRE") +
    ggtitle(paste("GLM ID link: 2017", enzTS17[j], sep = " "))
  ggsave(filename = paste("C:/Users/emily/OneDrive - University of Wyoming/Thesis Work/Thesis Writing/Figures/", "2017_GLMIDlink_", enzTS17[j], "_LinePlot.jpeg", sep = ""), plot = p,
         device = "jpeg")
}
