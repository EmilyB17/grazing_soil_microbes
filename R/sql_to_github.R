## Uploading data from SQL to GitHub
require(RODBC)
require(tidyverse)
dbConn<- odbcDriverConnect(connection = "driver={SQL Server Native Client 11.0};server=ECOINFO2,6003;database=Bean_SoilMicrobe;Uid=ebean;Pwd=Datadog1!")

## Enzymes
enz <- sqlQuery(dbConn, "SELECT        [2018luPlots].Plot, [2018luPlots].Treatment, [2018luPlots].Block, [2018luSamplePeriods].GrazeTime, [2018tblEnzymes_Vertical].Substrate, [2018tblEnzymes_Vertical].Enzyme_nm_g_hr 
FROM            [2018luLab_ID] INNER JOIN
                [2018luPlots] ON [2018luLab_ID].Plot = [2018luPlots].Plot INNER JOIN
                [2018luSamplePeriods] ON [2018luLab_ID].SampDate = [2018luSamplePeriods].SampDate INNER JOIN
                [2018tblEnzymes_Vertical] ON [2018luLab_ID].lab_id = [2018tblEnzymes_Vertical].lab_ID")
# write to GH
write.table(enz, file = "./data/2018enzymes_vertical.txt", sep = "\t", append = FALSE, row.names = FALSE)

# horizontal enzymes
enz.hor <- spread(data = enz, key = Substrate, value = Enzyme_nm_g_hr)
# write to GH
write.table(enz.hor, file = "./data/2018enzymes_horizontal.txt", sep = "\t", append = FALSE,
            row.names = FALSE)

# normalized enzymes for NMDS
# normalizing function
norm <- function(x){
  if(is.numeric(x) == FALSE){ 
    stop("Yo dummy, can only work numeric values")
  }
  # Let's calculate the normalized values
  return(x/max(x))
  
}

# loop to add normalized data 
for(i in 5:ncol(enz.hor)){
  # normalize current column values
  m<- norm(enz.hor[,i])
  # add new column
  enz.hor<- cbind(enz.hor, m)
  # rename the column
  names(enz.hor)[ncol(enz.hor)]<- names(enz.hor)[i]
}

# subset to have only normalized data
enz.norm <- enz.hor[, -(5:14)]
# write to GH
write.table(enz.norm, file = "./data/2018enzymes_RelAbundance.txt", sep = "\t", append = FALSE,
            row.names = FALSE)

## Environmental Data 
#### AS OF AUGUST 2018, HAVE MINERAL N & GRAVIMETRIC MOISTURE
env <- sqlQuery(dbConn, "SELECT        [2018luSamplePeriods].GrazeTime, [2018luPlots].Block, 
                         [2018luPlots].Treatment, [2018luPlots].Plot, [2018tblmineralN].no3_mgkgdrysoil, [2018tblmineralN].nh4_mgkgdrysoil, [2018tblmineralN].minN_mgkgdrysoil, [2018tblGravMoisture].grav_mois 
                FROM            [2018tblmineralN] INNER JOIN
                [2018luLab_ID] INNER JOIN
                [2018luPlots] ON [2018luLab_ID].Plot = [2018luPlots].Plot INNER JOIN
                [2018luSamplePeriods] ON [2018luLab_ID].SampDate = [2018luSamplePeriods].SampDate ON [2018tblmineralN].lab_ID = [2018luLab_ID].lab_id RIGHT OUTER JOIN
                [2018tblGravMoisture] ON [2018luLab_ID].lab_id = [2018tblGravMoisture].lab_ID")
# write to GH
write.table(env, file = "./data/2018biogeochemical_data.txt", sep = "\t", append = FALSE,
            row.names = FALSE)

## RPM
rpm <- sqlQuery(dbConn, "SELECT        [2018luLab_ID].Plot, [2018luPlots].Block, [2018luPlots].Treatment, [2018luSamplePeriods].GrazeTime, [2018tblRPM].reading_rpm, [2018tblRPM].biomass_kg_plot
FROM            [2018luLab_ID] INNER JOIN
                [2018luPlots] ON [2018luLab_ID].Plot = [2018luPlots].Plot INNER JOIN
                [2018luSamplePeriods] ON [2018luLab_ID].SampDate = [2018luSamplePeriods].SampDate INNER JOIN
                [2018tblRPM] ON [2018luLab_ID].SampDate = [2018tblRPM].SampDate AND [2018luLab_ID].Plot = [2018tblRPM].Plot")
write.table(rpm, file = "./data/2018rpm.txt", sep = "\t", append = FALSE, row.names = FALSE)

#RPM Trial
trial <- sqlQuery(dbConn, "SELECT        [2018tblrpmTrial].*
FROM            [2018tblrpmTrial]")
write.table(trial, file = "./data/2018rpm_fieldtrial.txt", sep = "\t", append = FALSE,
            row.names = FALSE)

## Mineral N, NPOC, DON
cn <- sqlQuery(dbConn, "SELECT        [2018luLab_ID].Plot, [2018luPlots].Block, [2018luPlots].Treatment, [2018luSamplePeriods].GrazeTime, [2018tblminN_NPOC_DON].NO3_mgkgdrysoil, [2018tblminN_NPOC_DON].NH4_mgkgdrysoil, 
                         [2018tblminN_NPOC_DON].NPOC_mgkgdrysoil, [2018tblminN_NPOC_DON].mineralN_mgkgdrysoil, [2018tblminN_NPOC_DON].DON_mgkgdrysoil
               FROM            [2018luLab_ID] INNER JOIN
               [2018luPlots] ON [2018luLab_ID].Plot = [2018luPlots].Plot INNER JOIN
               [2018luSamplePeriods] ON [2018luLab_ID].SampDate = [2018luSamplePeriods].SampDate RIGHT OUTER JOIN
               [2018tblminN_NPOC_DON] ON [2018luLab_ID].lab_id = [2018tblminN_NPOC_DON].Lab_ID")
write.table(cn, file = "./data/2018minN_NPOC_DON.txt", sep = "\t", append = TRUE)


