# Script to create time series definitions used for variable length experiments
# with ActivePower data; finding the baseline does not work for one machine, so we
# need to intervene manually
library(data.table)

source("PreprocessingUtility.R")

### Params ###

# Minmax normalization and machine selection do not play role in summary method
dataPath <- "../../data/2017_Jan_Aug_ActivePower1456VariableLengthMinmax/"
datasetFile <- "Wirkleistung_30_seconds_valid.rds"
minLength <- 1800
splitThreshold <- 1800

### Routine ###

# Standard procedure called in BatchUtility
timeSeriesOverviewTable <- summarizeNonBaseTimeSeriesForAggDataDir(inputDataDir = dataPath,
    attribute = "mean", splitThreshold = splitThreshold, minLength = minLength)

# Special treatment of AVT_04
# Display base values as table - AVT_04 has two one which differ a lot
# dataset <- readRDS(paste0(dataPath, "wirkleistung_30_seconds.rds"))
# dataset[, getBaseValues(mean), keyby = machine_name]
# Base for AVT_04 not satisfactory, manual intervention
avt04Data <- readRDS(paste0(dataPath, datasetFile))[machine_name == "AVT_04"]
machine04OverviewTable <- avt04Data[, summarizeNonBaseTimeSeriesForSingleSensor(oneSensorTable = .SD,
    attribute = "mean", splitThreshold = splitThreshold, minLength = minLength, sensorBaseValues = 0),
    by = .(machine_name, sensor_name)]
timeSeriesOverviewTable <- rbind(timeSeriesOverviewTable[machine_name != "AVT_04"],
                                 machine04OverviewTable)

# Remove very long time series (limit to one day)
timeSeriesOverviewTable <- timeSeriesOverviewTable[endTime_numeric - startTime_numeric <= 86400]

# Restore original order and save
setorder(timeSeriesOverviewTable, machine_name, startTime_numeric)
dir.create(paste0(dataPath, "results/"), showWarnings = FALSE)
saveRDS(timeSeriesOverviewTable, paste0(dataPath, "results/timeSeriesOverviewTableInput.rds"))
