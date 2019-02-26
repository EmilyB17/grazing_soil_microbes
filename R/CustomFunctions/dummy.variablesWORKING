


# create function to make dummy variables
dummy.variables <- function(data) {
  data$GrazePRE <- ifelse(data$GrazeTime %in% "PRE", 1, 0)
  data$Graze24H <- ifelse(data$GrazeTime %in% "24H", 1, 0)
  data$Graze1WK <- ifelse(data$GrazeTime %in% "1WK", 1, 0)
  data$Graze4WK <- ifelse(data$GrazeTime %in% "4WK", 1, 0)
  data$TrtHI <- ifelse(data$Treatment %in% "HI", 1, 0)
  data$TrtLO <- ifelse(data$Treatment %in% "LO", 1, 0)
  data$TrtNO <- ifelse(data$Treatment %in% "NO", 1, 0)
  data$Block1 <- ifelse(data$Block %in% 1, 1, 0)
  data$Block2 <- ifelse(data$Block %in% 2, 1, 0)
  data$Block3 <- ifelse(data$Block %in% 3, 1, 0)
  data$Block4 <- ifelse(data$Block %in% 4, 1, 0)
  data$GrazeTime <- NULL
  data$Treatment <- NULL
  data$Block <- NULL
  return(data)
}

df <- dummy.variables(cn)

# run GLM on the dummy variables model


# fit a Gamma distribution with a log link
(get.models(dredge(
  glm(grav_mois ~ .,
      family = Gamma(link = "log"), data = df, na.action = "na.fail")),
  subset = cumsum(weight) <= 0.95))[1]

# check the normality of the residuals with a QQplot
qqnorm(resid(glm(formula = grav_mois ~ GrazeTime + 1, family = Gamma(link = "log"), 
                 data = cn, na.action = "na.fail")))
qqline(resid(glm(formula = grav_mois ~ GrazeTime + 1, family = Gamma(link = "log"), 
                 data = cn, na.action = "na.fail")))

# summarize the best fit model
summary(glm(formula = grav_mois ~ GrazeTime + 1, family = Gamma(link = "log"), 
            data = cn, na.action = "na.fail"))