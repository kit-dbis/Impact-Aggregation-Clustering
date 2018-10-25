# Script to create time series definitions used for variable length experiments
# with PowerFactor data; finding the baseline does not work for some machines,
# so we need to intervene manually
library(data.table)

source("PreprocessingUtility.R")

### Params ###

dataPath <- "../../data/2017_Jan_Aug_PowerFactorWithout13VariableLength/"
datasetFile <- "Leistungsfaktor2_30_seconds_valid.rds"
minLength <- 1800
splitThreshold <- 1800

### Routine ###

dataset <- readRDS(paste0(dataPath, datasetFile))
# Multiple base values recognized for multiple machines; 1 is more reasonable from plots
# dataset[, getBaseValues(mean), keyby = machine_name]
timeSeriesOverviewTable <- dataset[, summarizeNonBaseTimeSeriesForSingleSensor(
  oneSensorTable = .SD, attribute = "mean", splitThreshold = splitThreshold,
  minLength = minLength, sensorBaseValues = 1),
  by = .(machine_name, sensor_name)]

# Remove very long time series (limit to one day)
timeSeriesOverviewTable <- timeSeriesOverviewTable[endTime_numeric - startTime_numeric <= 86400]

# Restore original order and save
setorder(timeSeriesOverviewTable, machine_name, startTime_numeric)
dir.create(paste0(dataPath, "results/"), showWarnings = FALSE)
saveRDS(timeSeriesOverviewTable, paste0(dataPath, "results/timeSeriesOverviewTableInput.rds"))
