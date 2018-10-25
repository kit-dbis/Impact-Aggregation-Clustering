# Script to remove days from MeanVoltage dataset during which the value is
# below a threshold. Sometimes voltage is 0 (sensor switched off) or has short,
# but big drops which distort clustering results.
library(data.table)
source("PreprocessingUtility.R")

dataPath <- "../../data/2017_Jan_Aug_MeanVoltage237/"

# Return the days where the "mean" attribute always is above threshold
findValidDaysMeanVoltage <- function(sensorDataTable) {
  return(sensorDataTable[, .(isValid = all(mean >= 390, na.rm = TRUE)),
                         by = day][isValid == TRUE, day])
}

preprocessAndReduceAggDataDir(dataPath = dataPath, validDayFunction = findValidDaysMeanVoltage)
