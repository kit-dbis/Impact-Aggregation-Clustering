# Script to remove days from PositiveEnergy dataset during which the value is
# decreasing. According to a domain expert, this sensor's values should only
# increase (in 1 kWh steps), but they fall to 0 (sensor switched off) or have
# short drops from time to time.
library(data.table)
source("PreprocessingUtility.R")

dataPath <- "../../data/2017_Jan_Aug_PositiveEnergy1234789Minmax/"

# Return the days where the "mean" attribute is monotonically increasing
findValidDaysEnergy <- function(sensorDataTable) {
  return(sensorDataTable[, .(isValid = all(mean >= shift(mean), na.rm = TRUE)),
                         by = day][isValid == TRUE, day])
}

preprocessAndReduceAggDataDir(dataPath = dataPath, validDayFunction = findValidDaysEnergy)
