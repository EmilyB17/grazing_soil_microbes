
## Year 2: Univariate Analyses

## ---- getData ----

## require packages
require(dplyr)
require(tidyr)
require(emmeans)
require(ggplot2)

# read prepData function
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/prepData2018.R")
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")
## read in data
# all C and N
cn <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018CN_data_updated10_18.txt", header = TRUE))

# enzymes
enz <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018enzymes_vertical.txt", header = TRUE))

# vegetation RPM and calculated biomass (not raw biomass - double check this)
veg <- prepData2018(read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/2018rpm.txt", header = TRUE))

# combine
all <- cn %>% 
  merge(veg, by = c("Plot", "GrazeTime", "Block", "Treatment"))

## REMOVE BLOCK 4 FROM ALL
all <- all %>% filter(!Block == 4)
veg <- veg %>% filter(!Block == 4)

## ---- explore ----

# make vertical for easy exploration
allv <- all %>% 
  gather(key = "Param", value = "value", grav_mois, DON_mgkgdrysoil, NPOC_mgkgdrysoil, mineralN_mgkgdrysoil, NH4_mgkgdrysoil, NO3_mgkgdrysoil,
         reading_rpm, biomass_kg_plot) %>% 
  mutate(Param = as.factor(Param))

# perform normality testing in a loop
params <- unique(allv$Param)
shapiros <- list() 

for(i in 1:length(params)) {
  
  # print histogram
  hist(allv$value[allv$Param == params[i]], main = paste0("histogram: ", params[i]))
  
  # print shapiro-wilkes test
  s <- shapiro.test(allv$value[allv$Param == params[i]])
  s$data.name <- as.character(params[i])
  shapiros <- c(shapiros, list(s))
  
  # print QQplot
  qqnorm(allv$value[allv$Param == params[i]], main = paste0("QQ plot: ", params[i]))
  qqline(allv$value[allv$Param == params[i]])
  
}

# view normality test
shapiros

# remove extreme NO3 outlier - we also lose this data point for mineral N (NO3 + NH4)
all$NO3_mgkgdrysoil[which(all$NO3_mgkgdrysoil > 80)] <- NA
all$mineralN_mgkgdrysoil[which(all$mineralN_mgkgdrysoil > 80)]  <- NA

# calculate percent differences and log transform
allpd<- percentDifferences(df = all,
                               ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                               timeKey = "GrazeTime",
                               timeLevels = c("PRE", "24H", "1WK", "4WK"),
                               level1 = "PRE")
allpdv <- allpd %>% 
  gather(key = "Param", value = "value", biomass_kg_plot, DON_mgkgdrysoil, grav_mois,
         mineralN_mgkgdrysoil, NH4_mgkgdrysoil, NO3_mgkgdrysoil, NPOC_mgkgdrysoil, reading_rpm) %>% 
  mutate(logvalue = log1p(value + 101)) %>% 
  select(-value)

# make horizontal version
allpdvh <- allpdv %>% 
  spread(key = Param, value = logvalue)

# look at normality
params <- unique(allpdv$Param)
shapiros <- list() 

for(i in 1:length(params)) {
  
  # print histogram
  hist(allpdv$logvalue[allpdv$Param == params[i]],  main = paste0("histogram: ", params[i]))
  
  # print shapiro-wilkes test
  s <- shapiro.test(allpdv$logvalue[allpdv$Param == params[i]])
  s$data.name <- as.character(params[i])
  shapiros <- c(shapiros, list(s))
  
  # print QQplot
  qqnorm(allpdv$logvalue[allpdv$Param == params[i]], main = paste0("QQ plot: ", params[i]))
  qqline(allpdv$logvalue[allpdv$Param == params[i]])
  
}

## ---- GLM CN data ----

allpdv$Param <- factor(allpdv$Param)
params <- unique(allpdv$Param)
outDF1 <- data.frame()
outDF2 <- data.frame()
modelFit <- data.frame()
for(i in 1:length(params)) {
  mod <- glm(logvalue ~ diffTimeSeries * Treatment,
             data = filter(allpdv, Param == params[i]),
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
    emmeans(mod, pairwise ~ Treatment | diffTimeSeries, type = "response")$contrasts)
  e$Param <- factor(as.character(params[i]))
  outDF1 <- rbind(e, outDF1)
  e2 <- as.data.frame(
    emmeans(mod, pairwise ~ diffTimeSeries | Treatment, type = "response")$contrasts)
  e2$Param <- factor(as.character(params[i]))
  outDF2 <- rbind(e2, outDF2)
  
}

outDF1$sign <- ifelse(outDF1$p.value < 0.05,
                      "significant",
                      "not_significant")

outDF2$sign <- ifelse(outDF2$p.value < 0.05,
                      "significant",
                      "not_significant")
sig.TS1 <- outDF1[which(outDF1$sign %in% "significant"),]
sig.Trt1 <- outDF2[which(outDF2$sign %in% "significant"),]

## ---- Significance Boxplots ----

sigparam <- unique(sig.TS1$Param)
sigparamTrt <- unique(sig.Trt1$Param)

### NPOC
# get significant contrasts
sigparam[2]
npoc <- sig.TS1[sig.TS1$Param == "NPOC_mgkgdrysoil",] # at 1WK, HI-LO; at 24H, HI-LO
npoc1 <- sig.Trt1[sig.Trt1$Param == "NPOC_mgkgdrysoil", ] # none

# make dataframe for line plots
datf <- allpdv %>% 
  group_by(Treatment, diffTimeSeries, Param) %>% 
  summarize(mean = mean(logvalue, na.rm = T),
            sd = sd(logvalue, na.rm = T),
            cilo = mean - 2*sd,
            cihi = mean + 2*sd) %>% 
  ungroup() %>% 
  gather(key = stat, value = value, mean, sd, cilo, cihi) %>% 
  spread(key = diffTimeSeries, value = value) %>% 
  mutate(PRE = 0) %>% 
  gather(key = diffTimeSeries, value = value, PRE, diff_24H, diff_1WK, diff_4WK) %>% 
  mutate(diffTimeSeries = factor(diffTimeSeries, ordered = TRUE, levels = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"))) %>% 
  spread(key = stat, value = value)

## DOES IT MATTER IF PRE IS CALCULATED IN THE MEANS OR NOT
datf1 <- allpdv %>% 
  spread(key = diffTimeSeries, value = logvalue) %>% 
  mutate(PRE = 0) %>% 
  gather(key = diffTimeSeries, value = logvalue, PRE, diff_24H, diff_1WK, diff_4WK) %>% 
  mutate(diffTimeSeries = factor(diffTimeSeries, ordered = TRUE, levels = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"))) %>% 
  group_by(Treatment, diffTimeSeries, Param) %>% 
  summarize(mean = mean(logvalue, na.rm = T),
            sd = sd(logvalue, na.rm = T),
            cilo = mean - 2*sd,
            cihi = mean + 2*sd) %>% 
  ungroup() 
  
  
  # contrasts: 
ggplot(data = datf[datf$Param %in% "NPOC_mgkgdrysoil", ], aes(x = diffTimeSeries, y = mean, group = Treatment,
                                                              color = Treatment)) +
  geom_point() +
  geom_line() 
