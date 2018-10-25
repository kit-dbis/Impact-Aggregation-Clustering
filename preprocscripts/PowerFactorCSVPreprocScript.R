# Script to remove days from PowerFactor dataset during which the value is smaller
# than 1 (base level) for less than half an hour or where the value is smaller
# than 1 at beginning or end of the time series (no on/off behaviour during day).
# Routine quite similar to MeanAmperage pre-processing (which has a different base).
library(data.table)
source("PreprocessingUtility.R")

dataPath <- "../../data/2017_Jan_Aug_PowerFactorWithout13/"

# Return the days where the "mean" attribute is 1 for time series begin and
# end (base value, machines are switched on/off during day), but less than 1
# for at least half an hour (remove time series with outlier spikes)
findValidDaysPowerFactor <- function(sensorDataTable) {
  return(sensorDataTable[, .(isValid = sum(mean < 1) >= 100 & mean[1] == 1 &
                               mean[.N] == 1), by = day][isValid == TRUE, day])
}

preprocessAndReduceAggDataDir(dataPath = dataPath, validDayFunction = findValidDaysPowerFactor)
