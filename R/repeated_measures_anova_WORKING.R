# ---- Repeated Measures ANOVA ----
## Emily Bean

# Between-subjects variable: Treatment (each plot has a different treatment)
# Within-subjects variable: GrazeTime (each plot repeated four GrazeTimes)

require(ggplot2)
require(tidyverse)
require(dplyr)
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/prepData.R")

# ---- Microbial Biomass C ----
toc <- read.table(file = "toc.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
mbc <- toc[, -(6:20)]
mbc <- mbc[, -(7:8)]

# Examine data relative to the PRE baseline
pre <- aggregate(mbc_mgg ~ Plot, data = mbc[mbc$GrazeTime %in% "PRE",], mean)
pre$mean <- pre$mbc_mgg
pre$mbc_mgg <- NULL
dat <- aggregate(mbc_mgg ~ GrazeTime + Plot + Treatment, data = mbc[!mbc$GrazeTime %in% "PRE",], mean) 

all <- merge(pre, dat, by = "Plot")
all$delta <- all$mean - all$mbc_mgg
mbc <- all


# Now the data is averaged by Plot, and "delta" is the difference of the average of that plot from the average of the same plot at the PRE baseline

# plot
ggplot(data = all, aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# perform repeated-measures ANOVA
mod.mbc <- aov(all$delta ~ all$Treatment * all$GrazeTime + Error(all$delta / all$GrazeTime))
summary(mod.mbc)

# Holm procedure
with(mbc, pairwise.t.test(delta, Treatment, paired = T)) # no significance
with(mbc, pairwise.t.test(delta, GrazeTime, paired = T)) #24H and 4WK are significant

# Bonferroni test - is the most conservative
with(mbc, pairwise.t.test(delta, Treatment, p.adjust.method = "bonferroni", paired = TRUE))
with(mbc, pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # 24H and 4Wk are significant

# Tukey's HSD test
(TukeyHSD(mod.mbc, conf.level = 0.95))

# ---- Microbial Biomass N ----
toc <- read.table(file = "toc.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
mbn <- toc[, -(6:21)]
mbn <- mbn[, -(7)]




# Examine data relative to the PRE baseline
pre <- aggregate(mbn_mgg ~ Plot, data = mbn[mbn$GrazeTime %in% "PRE",], mean)
pre$mean <- pre$mbn_mgg
pre$mbn_mgg <- NULL
dat <- aggregate(mbn_mgg ~ GrazeTime + Plot + Treatment, data = mbn[!mbn$GrazeTime %in% "PRE",], mean) 

all <- merge(pre, dat, by = "Plot")
all$delta <- all$mean - all$mbn_mgg
mbn <- all


# Now the data is averaged by Plot, and "delta" is the difference of the average of that plot from the average of the same plot at the PRE baseline

# plot
ggplot(data = all, aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# perform repeated-measures ANOVA
mod.mbn <- aov(all$delta ~ all$Treatment * all$GrazeTime + Error(all$delta / all$GrazeTime))
summary(mod.mbn) # GrazeTime is significant

# Holm procedure
with(mbn, pairwise.t.test(delta, Treatment, paired = T)) # no significance
with(mbn, pairwise.t.test(delta, GrazeTime, paired = T)) #24H-1WK, 24H-4WK, 1Wk-4WK are significant

# Bonferroni test - is the most conservative
with(mbn, pairwise.t.test(delta, Treatment, p.adjust.method = "bonferroni", paired = TRUE)) # no significance
with(mbn, pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # All are significant

# ---- NPOC ----
toc <- read.table(file = "toc.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
npoc <- toc[, -(6:12)]
npoc <- npoc[, -(7:16)]

# Examine data relative to the PRE baseline
pre <- aggregate(NPOC_t0Nmin_mgg ~ Plot, data = npoc[npoc$GrazeTime %in% "PRE",], mean)
pre$mean <- pre$NPOC_t0Nmin_mgg
pre$NPOC_t0Nmin_mgg <- NULL
dat <- aggregate(NPOC_t0Nmin_mgg ~ GrazeTime + Plot + Treatment, data = npoc[!npoc$GrazeTime %in% "PRE",], mean) 

all <- merge(pre, dat, by = "Plot")
all$delta <- all$mean - all$NPOC_t0Nmin_mgg
npoc <- all


# Now the data is averaged by Plot, and "delta" is the difference of the average of that plot from the average of the same plot at the PRE baseline

# plot
ggplot(data = all, aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# perform repeated-measures ANOVA
mod.npoc <- aov(all$delta ~ all$Treatment * all$GrazeTime + Error(all$delta / all$GrazeTime))
summary(mod.npoc)

# Holm procedure
with(npoc, pairwise.t.test(delta, Treatment, paired = T)) # no significance
with(npoc, pairwise.t.test(delta, GrazeTime, paired = T)) # no significance

# Bonferroni test - is the most conservative
with(npoc, pairwise.t.test(delta, Treatment, p.adjust.method = "bonferroni", paired = TRUE)) # no significance
with(npoc, pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # no significance

# ---- DON ----
toc <- read.table(file = "toc.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
don <- toc[, -(6:13)]
don <- don[, -(7:15)]

# Examine data relative to the PRE baseline
pre <- aggregate(DON_t0Nmin_mgg ~ Plot, data = don[don$GrazeTime %in% "PRE",], mean)
pre$mean <- pre$DON_t0Nmin_mgg
pre$DON_t0Nmin_mgg <- NULL
dat <- aggregate(DON_t0Nmin_mgg ~ GrazeTime + Plot + Treatment, data = don[!don$GrazeTime %in% "PRE",], mean) 

all <- merge(pre, dat, by = "Plot")
all$delta <- all$mean - all$DON_t0Nmin_mgg
don <- all


# Now the data is averaged by Plot, and "delta" is the difference of the average of that plot from the average of the same plot at the PRE baseline

# plot
ggplot(data = all, aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# perform repeated-measures ANOVA
mod.don <- aov(all$delta ~ all$Treatment * all$GrazeTime + Error(all$delta / all$GrazeTime))
summary(mod.don)

# Holm procedure
with(don, pairwise.t.test(delta, Treatment, paired = T)) # no significance
with(don, pairwise.t.test(delta, GrazeTime, paired = T)) #1WK-4WK and 24H-4WK are significant

# Bonferroni test - is the most conservative
with(don, pairwise.t.test(delta, Treatment, p.adjust.method = "bonferroni", paired = TRUE)) # no significance
with(don, pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # 1WK-4Wk and 24H-4Wk are significant

# ---- Enzymes ----
# Read in data
enz.vert <- prepData(read.table(file = "https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/enzymes_vertical.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE))


# Examine data relative to the PRE baseline
pre <- aggregate(Enzyme_nm_g_hr ~ Plot + Substrate, data = enz.vert[enz.vert$GrazeTime %in% "PRE",], mean)
pre$pre_mean <- pre$Enzyme_nm_g_hr
pre$Enzyme_nm_g_hr <- NULL
dat <- aggregate(Enzyme_nm_g_hr ~ GrazeTime + Plot + Treatment + Substrate, 
                 data = enz.vert[!enz.vert$GrazeTime %in% "PRE",], mean) 

all <- merge(pre, dat, by = c("Plot", "Substrate"))
all$delta <- all$Enzyme_nm_g_hr - all$pre_mean
enz.vert <- all
enz.vert$pre_mean <- NULL
enz.vert$Enzyme_nm_g_hr <- NULL

# Make horizontal
#enz <- spread(data = enz.vert, key = Substrate, value = delta)

# Now the data is averaged by Plot, and "delta" is the difference of the average of that plot from the average of the same plot at the PRE baseline


# Loop to look at each substrate in the repeated-measures ANOVA format
sub <- unique(enz.vert$Substrate)
# save the model objects
mods<- list()
for(i in 1:length(sub)) {
  mods<- c(mods, list(aov(delta ~ Treatment * GrazeTime + Error(delta / GrazeTime), 
                          data = filter(enz.vert, Substrate == sub[i]))))

}
names(mods)<- sub
lapply(mods, summary)

# CBH, NAG, PER1, PHOS: GrazeTime
# Do individual post-hoc tests for the substrates that showed significance with GrazeTime

## CBH
# plot
ggplot(data = enz.vert[enz.vert$Substrate %in% "CBH",], 
       aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# Holm procedure
with(enz.vert[enz.vert$Substrate %in% "CBH",], 
     pairwise.t.test(delta, GrazeTime, paired = T)) # 24H-1WK, 24H-4WK

# Bonferroni test - is the most conservative
with(enz.vert[enz.vert$Substrate %in% "CBH",], 
     pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) #24H-1WK, 24H-4WK

## NAG
# plot
ggplot(data = enz.vert[enz.vert$Substrate %in% "NAG",], 
       aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# Holm procedure
with(enz.vert[enz.vert$Substrate %in% "NAG",], 
     pairwise.t.test(delta, GrazeTime, paired = T)) # No significance, interesting since aov showed significance with GrazeTime

# Bonferroni test - is the most conservative
with(enz.vert[enz.vert$Substrate %in% "NAG",], 
     pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # Still no significance

## PER1
# plot
ggplot(data = enz.vert[enz.vert$Substrate %in% "PER1",], 
       aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# Holm procedure
with(enz.vert[enz.vert$Substrate %in% "PER1",], 
     pairwise.t.test(delta, GrazeTime, paired = T)) # No significance

# Bonferroni test - is the most conservative
with(enz.vert[enz.vert$Substrate %in% "PER1",], 
     pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # Still no significance

## PHOS
# plot
ggplot(data = enz.vert[enz.vert$Substrate %in% "PHOS",], 
       aes(x = Treatment, y = delta, colour = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ GrazeTime)

# Holm procedure
with(enz.vert[enz.vert$Substrate %in% "PHOS",], 
     pairwise.t.test(delta, GrazeTime, paired = T)) # 24H-1WK, 1WK-4WK

# Bonferroni test - is the most conservative
with(enz.vert[enz.vert$Substrate %in% "PHOS",], 
     pairwise.t.test(delta, GrazeTime, p.adjust.method = "bonferroni", paired = TRUE)) # 24H-1WK, 1WK-4WK



