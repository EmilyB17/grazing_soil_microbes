# Extracellular enzymes NMDS script
# Creates ordination in vegan, then plots in ggplot

# //// as tutorial, only run all commands (don't change anything)
# to help your understanding, check the structure of variables in the console as you create them
# example: 
x <- 1
str(x)
# //// BEFORE MANIPULATING THIS SCRIPT IN R:
# Save to your own directory
# set your own working directory

setwd()

# part of this tutorial comes from: https://oliviarata.wordpress.com/2014/04/17/ordinations-in-ggplot2/
# part of this tutorial comes from Van Diepen UW lab
require(vegan)
require(ggplot2)
require(grid) # for the arrows in ggplot, I have not looked into this package further

# //// You need 2 dataframes: grouping & environmental variables, and normalized enzyme data 

# environmental data - this includes grouping variables as well as environmental variables
# PERMANOVA (adonis function) can only examine factors - make sure all variables are factors
# my grouping variables are Plot, Sample, Block, Treatment, and GrazeTime 
# my environmental variables are ammonium.ppm and gravimetric.moisture
dat <- read.table(file = "https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/environmental_data.txt",
                  sep = "\t", header = TRUE)
dat$ammonium.ppm <- factor(dat$ammonium.ppm)
dat$gravimetric.moisture <- factor(dat$gravimetric.moisture)
dat$Block <- factor(dat$Block)
dat$bd <- NULL
dat$ph <- NULL
dat$Plot <- NULL
dat$Sample <- NULL


# enzyme data (aka: abundance data or species data)
# this data is already normalized to relative abundance
allgrz <- read.table(file = "https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/enzymes_relativeabundance.txt",
                     sep = "\t", header = TRUE)
allgrz$Plot <- NULL
allgrz$Block <- NULL
allgrz$Treatment <- NULL
allgrz$Sample <- NULL
allgrz$GrazeTime <- NULL
str(allgrz)

# function metaMDS() to perform an ordination on relative abundance data
ord <- metaMDS(allgrz)

# function vegdist() to calculate a Bray-Curtis distance matrix of relative abundance data
dis <- vegdist(allgrz)

# calculate a stressplot to examine how well your ordination fits the actual dissimilarity data
# if the stressplot R2 is low, this means your dissimilarity data does not fit well into 
# a 2 dimensional ordination. 
stressplot(ord, dis)

# adonis() function calculates a permanova on the environmental and/or grouping variables
adonis(dis ~ ammonium.ppm + gravimetric.moisture + trt.grz,
       data = dat, permutations = 999)

# the envfit() function fits your environmental and/or grouping data onto the ordination
# this gives p values that you can later use to plot onto the ordination 
# as correlating environmental variables
eft <- envfit(ord ~., data = dat, perm = 1000, na.rm = TRUE)
eft

# This is the end of creating a basic ordination
# Quickly we will plot this in vegan
plot(ord, type = "p", display = "sites") # all of the ordination points
plot(ord, type = "text", display = "species", lty = 1) # plots enzymes 
with(dat, ordiellipse(ord, Treatment, kind = "se", conf = 0.95, label = TRUE)) # adds ellipses by grouping variable 
# with 95% confidence intervals
plot(eft, p.max=0.05, col="black", cex=0.8) # plots significant environmental variables from the envfit() function

# Now, we begin prepping to plot the ordination in ggplot
# since ggplot cannot plot ordinations, we need to pull out each individual part

# First pull the coordinates from the metaMDS() ordination into 
# your grouping & environmental variables dataframe
dat$NMDS1 <- ord$points[,1]
dat$NMDS2 <- ord$points[,2]
head(dat) # you should see NMDS1 and NMDS2 as columns 
# this allows us to plot all rows of data in ordination form

# This pulls x and y coordinates for each 
# abundance variable - in this case, enzymes (other cases, taxa or species)
species.scores <- data.frame(scores(ord, display = "species"))
species.scores$species <- row.names(species.scores) # if you have rownames this imports them
# this will be used later to plot the enzymes on their ordination coordinates

# function for ellipses that is floating all around the internet
# only run this function, don't change anything in it
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# data for ellipses
# loop: for each level of your grouping variable, save the correct coordinates for an ellipse 
# in the dataframe df_ell
# my grouping variable is Treatment and my environmental variables dataframe is dat

df_ell <- data.frame() # creates empty dataframe to fill in the loop
for(g in levels(dat$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dat [dat$Treatment==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,Treatment=g))
}
str(df_ell) # your dataframe should have 3 columns: NMDS1, NMDS2, and your grouping variable

# This pulls data to label the centers of the ellipses
# dat is the environmental/grouping variables dataframe
NMDS.mean = aggregate(dat[,c("NMDS1", "NMDS2")],
                      list(group = dat$Treatment), mean)


# plot in ggplot
# if it bothers you to have it all under one command (GORDON), 
# add a variable in front of ggplot() and in front of every line 
# or comment out lines one at a time to work on them

# hint: everywhere you see "Treatment", you will need to replace with your own grouping variable
ggplot(data = dat, aes(x = NMDS1, y = NMDS2)) + # this sets the basic plot
  geom_polygon(data = df_ell, aes(fill = Treatment, group = Treatment, alpha=Treatment)) + #add ellipses, grouped by your grouping variable
  # if you don't want the polygon filled with a color, use geom_path instead of geom_polygon
  geom_point(aes(x = NMDS1, y = NMDS2, shape = Treatment, color = Treatment), size = 1) +  #adds ordination points, shaped by grouping variable
  annotate("text",x = NMDS.mean$NMDS1,y = NMDS.mean$NMDS2,label=NMDS.mean$group) + #labels for the centroids of the ellipses
  geom_segment(data = species.scores,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), # creates segment for enzymes
               arrow = arrow(length = unit(0.20, "cm")), colour = "grey") +  # creates arrow on the end of the segment
  geom_text(data = species.scores, # labels the enzymes at the end of the arrows 
            aes(x = NMDS1, y = NMDS2, label=species),
            size = 3)  +
  theme_classic() + # white background, no grid lines, black axises
  # theme_bw() will keep gridlines but put thick black line around plot
  ggtitle("Enzyme Relative Abundance NMDS, grouped by Treatment")
