# Script to remove days from MeanAmperage dataset during which the value is greater
# than 0 (base level) for less than half an hour or where the value is greater than
# 0 at beginning or end of the time series (no on/off behaviour during day).
library(data.table)
source("PreprocessingUtility.R")

dataPath <- "../../data/2017_Jan_Aug_MeanAmperageWithout3Minmax/"

# Return the days where the "mean" attribute is zero for time series begin and
# end (switched on/off during day), but greater than zero for at least half an
# hour (remove time series with outlier spikes)
findValidDaysMeanAmperage <- function(sensorDataTable) {
  return(sensorDataTable[, .(isValid = sum(mean > 0) >= 100 & mean[1] == 0 &
      mean[.N] == 0), by = day][isValid == TRUE, day])
}

preprocessAndReduceAggDataDir(dataPath = dataPath, validDayFunction = findValidDaysMeanAmperage)
