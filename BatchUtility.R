# Contains functions for batch scripts

source("AggregationUtility.R")
source("ClusteringEvaluationUtility.R")
source("ClusteringUtility.R")
source("DissimilarityUtility.R")
source("PreprocessingUtility.R")

library(batchtools)
library(data.table)
library(optparse)
library(microbenchmark) # add to Packrat
library(parallel)

# Calls a function which works with a batchtools registry (parameter "expRegistry"
# necessary), measures time and cleans up if successful. All parameters of the
# function (apart from the registry) should be passed as named parameters.
callBatchtoolsRoutine <- function(jobName, FUN, ...) {
  cat("---Starting batchtools routine \"", jobName, "\"---\n", sep = "")
  experimentRegistry <- makeExperimentRegistry(file.dir = "registry", seed = 1)
  params <- list(...)
  params[["expRegistry"]] <- experimentRegistry
  if (is.null(params[["cpuCores"]])) {
    params[["cpuCores"]] <- detectCores()
  }
  startTime <- as.numeric(Sys.time())
  result <- do.call(FUN, params)
  waitForJobs(reg = experimentRegistry)
  endTime <- as.numeric(Sys.time())
  cat("The routine \"", jobName, "\" took ", round(endTime - startTime), "s\n", sep = "")
  if (nrow(findErrors(reg = experimentRegistry)) != 0) {
    saveRDS(result, "batchtoolsResultIncomplete.rds")
    stop(paste0("Errors occured during \"", jobName, "\"."))
  } else {
    removeRegistry(reg = experimentRegistry, wait = 0)
  }
  gc()
  return(result)
}

# Returns a list of dissimilarities, considering if aggregations are uni- or
# multivariate and if time series lengths are fixed or variable
getDefaultDissimilarities <- function(aggNamesList, intervalLength) {
  result <- if (all(sapply(aggNamesList, length) == 1)) {
    DEFAULT_DISSIMILARITIES
  } else {# raw aggregation at least univariate, at least one multivariate
    DEFAULT_DISPATCH_DISSIMILARITIES
  }
  if (intervalLength <= 0) { # only variable length dissimilarities
    # Lp norms + CORT assume equal length time series; PDD is undefined if one
    # ts length == 1 (which might happen if short time series aggregated);
    # DTW.Band would need proper windowSize parameter (there might be no solution
    # if big length difference)
    result <- result[intersect(names(result),
                               c("DTW", "DTW.CID", "EDR", "ERP", "SBD", "CDM"))]
    for (disName in c("DTW", "DTW.CID", "EDR", "ERP")) {
      if (disName %in% names(result)) {
        result[[disName]]$params <- c(result[[disName]]$params, normalize = TRUE)
      }
    }
  }
  return(result)
}

# Returns common time series lengths in energy data research, considering one day
# is made up of intervals with the "intervalLength" (in seconds, default: each day
# represented by only one interval)
getDefaultDailyAggregationLevels <- function(intervalLength = 86400) {
  if (intervalLength <= 0) { # negative lengths represent variable length time
    # series, we use default levels (which will be be understood in a relative
    # rather than absolute way during aggregation)
    return(c(1440, 288, 144, 96, 48, 24, 12, 4))
  }
  intervalsPerDay <- 86400 / intervalLength
  if (intervalsPerDay != floor(intervalsPerDay)) {
    stop("Your interval length does not divide the length of one day without remainder.")
  }
  # 1 min, 5 min, 10 min, 15 min, 30 min, 1 h, 2 h, 6 h 
  result <- c(1440, 288, 144, 96, 48, 24, 12, 4) / intervalsPerDay
  result <- result[result >= 2] # keep only levels with at least 2 points
  if (any(result != floor(result))) {
    stop("Your interval length is not compatible with the default aggregation levels.")
  }
  return(result)
}

# Takes a string with the aggregation names and creates a list
# - Substrings separated by pipes (|) become elements of the list (aggregation types).
# - Elements might be multivariate if separated by an underscore (_), they are
# converted to character vectors.
# Whitespaces are removed.
convertAggregationStringToList <- function(aggregationString) {
  aggregationString <- gsub("[[:space:]]", replacement = "", x = aggregationString)
  result <- strsplit(aggregationString, split = "|", fixed = TRUE)[[1]]
  result <- strsplit(result, split = "_", fixed = TRUE)
  return(result)
}

# Takes a string with machine names and creates a character vector
# Substrings separated by pipes (|) become elements of the vector.
# Whitespaces are removed.
convertMachineStringToVector <- function(machineString) {
  machineString <- gsub("[[:space:]]", replacement = "", x = machineString)
  result <- strsplit(machineString, split = "|", fixed = TRUE)[[1]]
  return(result)
}

# Checks if there is only one *.rds file (representing aggregated datasets) under
# the "aggDataPath" and if yes, uses the functions in "aggFunctionsList" to create
# aggregates locally
saveLocalAggregatesIfNecessary <- function(datasetName, aggDataPath, aggFunctionsList) {
  aggDataTimeSeriesFiles <- list.files(aggDataPath, pattern = "\\.rds", full.names = TRUE)
  if (length(aggDataTimeSeriesFiles) == 1) {
    rawDatasetList <- list(readRDS(aggDataTimeSeriesFiles[1]))
    names(rawDatasetList) <- datasetName
    outputBasePathList <- list(aggDataPath)
    names(outputBasePathList) <- datasetName
    callBatchtoolsRoutine("aggregation", aggregateAndSaveData, tsListList = rawDatasetList,
        aggList = aggFunctionsList, outputBasePathList = outputBasePathList)
  }
  return(TRUE)
}

# Saves already pre-processed aggregated datasets from Spark as time series
# lists, using only machines, sensors and days contained in the data.table
# "machineSensorDayCombinations"
# params: see "savePreprocessedTimeSeries()"
saveTimeSeriesForCombinations <- function(datasetName, dataPath,
    aggDataPath, machineSensorDayCombinations, intervalLength = 86400,
    aggNamesList = NULL, aggFunctionsList = NULL, normalize = "") {
  # Convert time series and save them, missing values are interpolated
  saveAllAggregateTimeSeries(datasetName = datasetName, inputDataDir = dataPath,
      outputDataDir = aggDataPath, machineDaySensorCombinations = machineSensorDayCombinations,
      intervalLength = intervalLength, aggColumnList = aggNamesList, rawAggColumn = "mean",
      naReplace = 0, normalize = normalize)
  # If only one database aggregation level supplied: create further levels with batchtools
  saveLocalAggregatesIfNecessary(datasetName = datasetName, aggDataPath = aggDataPath,
                                 aggFunctionsList = aggFunctionsList)
}

# Saves time series files for all aggregations and levels, using time intervals
# of fixed length (one day or intervals of a day). Before coverting table datasets
# to time series lists, it is checked which time intervals (per machine and sensor)
# are complete to a certain extend (number of timestamps) and which only contain
# the same value (i.e. machine/sensor off). Constant series and series with a certain
# amount of missing timestamps won't be created.
# parameters: see "savePreprocessedTimeSeries()"
saveFixedLengthTimeSeries <- function(datasetName, dataPath, aggDataPath, resultsPath,
    machineNames = NULL, intervalLength = 86400, aggNamesList = NULL, aggFunctionsList = NULL,
    normalize = "", testAggDataPath = NULL, fullForecastAggDataPath = NULL,
    fullForecastTestAggDataPath = NULL, testRatio = 0) {
  # Determine "high quality" days - non-constant time series, only a few missing timestamps
  machineSensorDayCountTable <- summarizeValueCountForAggDataDir(dataPath, intervalLength = intervalLength)
  saveRDS(machineSensorDayCountTable, paste0(resultsPath, "machineSensorDayCountTable.rds"))
  machineSensorDayCombinations <- machineSensorDayCountTable[NRatio > 0.99,
       .(machine_name, sensor_name, day, timeOfDay)]
  machineSensorDayConstantTable <- summarizeConstantSeriesForAggDataDir(dataPath,
       attribute = "mean", intervalLength = intervalLength)
  machineSensorDayCombinations <- merge(machineSensorDayCombinations,
      machineSensorDayConstantTable[constant == FALSE,
                                    .(machine_name, sensor_name, day, timeOfDay)])
  if (!is.null(machineNames) & !all(machineNames == "")) {
    machineSensorDayCombinations <- machineSensorDayCombinations[machine_name %in% machineNames]
  }
  if (testRatio == 0) {# no forecasting evaluation, one output directory
    saveRDS(machineSensorDayCombinations, paste0(resultsPath, "machineSensorDayCombinations.rds"))
    saveTimeSeriesForCombinations(datasetName = datasetName, dataPath = dataPath,
        aggDataPath = aggDataPath, machineSensorDayCombinations = machineSensorDayCombinations,
        intervalLength = intervalLength, aggNamesList = aggNamesList,
        aggFunctionsList = aggFunctionsList, normalize = normalize)
  } else {# Also create a test set for forecasting evaluation and split input+target
    machineSensorDayCombinations <- removeIntervalsWithoutPredAndSucc(
      machineSensorDayCombinations, intervalLength = intervalLength)
    lastTrainDay <- machineSensorDayCombinations[, sort(day)][
      round((1 - testRatio) * nrow(machineSensorDayCombinations))]
    trainDataCombinations <- machineSensorDayCombinations[day <= lastTrainDay]
    testDataCombinations <- machineSensorDayCombinations[day > lastTrainDay]
    trainDataCombinations <- removeIntervalsWithoutPredAndSucc(trainDataCombinations,
        intervalLength = intervalLength)
    testDataCombinations <- removeIntervalsWithoutPredAndSucc(testDataCombinations,
        intervalLength = intervalLength)
    orderForAggDataTimeSeriesListConversion(trainDataCombinations)
    orderForAggDataTimeSeriesListConversion(testDataCombinations)
    saveRDS(trainDataCombinations, paste0(resultsPath, "trainMachineSensorDayCombinations.rds"))
    saveRDS(testDataCombinations, paste0(resultsPath, "testMachineSensorDayCombinations.rds"))
    # Training: input + target series (which are necessary for forecasting)
    cat("-Training set, full (input + target)-\n")
    saveTimeSeriesForCombinations(datasetName = datasetName, dataPath = dataPath,
        aggDataPath = fullForecastAggDataPath,
        machineSensorDayCombinations = trainDataCombinations,
        intervalLength = intervalLength, aggNamesList = aggNamesList,
        aggFunctionsList = aggFunctionsList, normalize = normalize)
    # Test: input + target series (which are necessary for forecasting)
    cat("-Test set, full (input + target)-\n")
    saveTimeSeriesForCombinations(datasetName = datasetName, dataPath = dataPath,
        aggDataPath = fullForecastTestAggDataPath,
        machineSensorDayCombinations = testDataCombinations,
        intervalLength = intervalLength, aggNamesList = aggNamesList,
        aggFunctionsList = aggFunctionsList, normalize = normalize)
    # Predecessor-successor relationship in training and test set (for forecasting)
    trainDataPrevActList <- list(createPrevActListForAggDataSummary(trainDataCombinations,
        intervalLength = intervalLength))
    names(trainDataPrevActList) <- datasetName
    testDataPrevActList <- list(createPrevActListForAggDataSummary(testDataCombinations,
        intervalLength = intervalLength))
    names(testDataPrevActList) <- datasetName
    saveRDS(trainDataPrevActList, paste0(resultsPath, "trainDataPrevActList.rds"))
    saveRDS(testDataPrevActList, paste0(resultsPath, "testDataPrevActList.rds"))
    # Limit to input series of forecasting (only these need to be clustered/assigned)
    trainDataCombinations <- trainDataCombinations[trainDataPrevActList[[datasetName]]$prev]
    testDataCombinations <- testDataCombinations[testDataPrevActList[[datasetName]]$prev]
    # Training: input series (which are clustered and used for assignment)
    cat("-Training set, input only-\n")
    saveTimeSeriesForCombinations(datasetName = datasetName, dataPath = dataPath,
        aggDataPath = aggDataPath, machineSensorDayCombinations = trainDataCombinations,
        intervalLength = intervalLength, aggNamesList = aggNamesList,
        aggFunctionsList = aggFunctionsList, normalize = normalize)
    # Test: input series (which are used for assignment)
    cat("-Test set, input only-\n")
    saveTimeSeriesForCombinations(datasetName = datasetName, dataPath = dataPath,
        aggDataPath = testAggDataPath, machineSensorDayCombinations = testDataCombinations,
        intervalLength = intervalLength, aggNamesList = aggNamesList,
        aggFunctionsList = aggFunctionsList, normalize = normalize)
  }
  return(TRUE)
}

# Saves time series files for all aggregations and levels, using variable lengths
# (split between time series determined by "intervalLength").
# parameters: see "savePreprocessedTimeSeries()"
saveVariableLengthTimeSeries <- function(datasetName, dataPath, aggDataPath, resultsPath,
    machineNames = NULL, intervalLength = 300, aggNamesList = NULL,
    aggFunctionsList = NULL, normalize = "") {
  if (intervalLength == 0) { # use pre-defined table with time series starts and ends
    inputOverviewTablePath <- paste0(resultsPath, "timeSeriesOverviewTableInput.rds")
    if (!file.exists(inputOverviewTablePath)) {
      stop(paste0("Using an intervalLength of 0 requires a file 'timeSeriesOverviewTableInput.rds' ",
                  "placed at the 'resultsPath'."))
    }
    timeSeriesOverviewTableInput <- readRDS(inputOverviewTablePath)
    # Only problem: overview table might have different sensor names, so we make
    # sure sensors are unique per machine and join based on machines
    machineSensorCombis <- summarizeMachinesAndSensorsForAggDataDir(inputDataDir = dataPath)
    if (any(machineSensorCombis[, uniqueN(sensor_name), by = machine_name]$V1 != 1)) {
      stop(paste0("Your dataset contains more than one sensor per machine. This ",
                  "is not supported in combination with 'intervalLength' == 0."))
    }
    if (any(timeSeriesOverviewTableInput[, uniqueN(sensor_name), by = machine_name]$V1 != 1)) {
      stop(paste0("Your file 'timeSeriesOverviewTableInput.rds' contains more ",
                  "than one sensor per machine. This is not supported in combination ",
                  "with 'intervalLength' == 0."))
    } else {
      timeSeriesOverviewTableInput[, sensor_name := NULL]
    }
    timeSeriesOverviewTable <- merge(machineSensorCombis, timeSeriesOverviewTableInput,
                                     by = "machine_name")
  } else {# interpret intervalLength as minimum "pause" between time series
    splitInterval <- max(1, abs(intervalLength))
    timeSeriesOverviewTable <- summarizeNonBaseTimeSeriesForAggDataDir(inputDataDir = dataPath,
        attribute = "mean", splitThreshold = splitInterval, minLength = 1800)
  }
  if (!is.null(machineNames) & !all(machineNames == "")) {
    timeSeriesOverviewTable <- timeSeriesOverviewTable[machine_name %in% machineNames]
  }
  saveRDS(timeSeriesOverviewTable, paste0(resultsPath, "timeSeriesOverviewTable.rds"))
  saveAllAggregateTimeSeries(datasetName = datasetName, inputDataDir = dataPath,
      outputDataDir = aggDataPath, timeSeriesCombinations = timeSeriesOverviewTable,
      intervalLength = intervalLength, aggColumnList = aggNamesList, rawAggColumn = "mean",
      naReplace = 0, normalize = normalize)
  # If only one database aggregation level supplied: create further levels with batchtools
  saveLocalAggregatesIfNecessary(datasetName = datasetName, aggDataPath = aggDataPath,
                                 aggFunctionsList = aggFunctionsList)
}

# Preprocesses CSV files produced by Spark aggregation jobs. CSVs are turned into RDS
# (with some column renaming and datatype conversion) and finally time series representing
# single days of machines and sensors, retaining only series which contain enough timestamps
# and are not constant.
# - "datasetName": used for the time series list *.rds files; must not contain underscores if
# used for further analysis
# - "dataPath": directory containing the CSV files from Spark; if multiple ones, they are
# used for creating the aggregate data sets; if only one, aggregates are computed locally
# - "aggDataPath": directory to put the time series list files used for clustering
# - "resultsPath": directory to put two info files which summarize machines, sensors
# and days contained in the datasets
# - "machineNames": character vector denoting machines to be used; if NULL (default)
# or empty character, all machines from the input dataset are considered
# - "intervalLength": if posititve: length (in seconds) of time intervals represented
# by one time series, default is one day (86400), but you can also use a divisor
# of that; negative: duration used to split time series (if at least "intervalLength"
# seconds the sensor base value (e.g. 0) occurs); if zero, a *.rds file is read
# which determines start and end times (also variable length)
# - "aggNamesList": list of character vectors of column names which will be
# extracted from aggegate tables (can be NULL if local aggregation should be
# performed)
# - "aggFunctionList": list of aggregation functions for local aggregation (can
# be NULL if aggregated datasets used as input)
# - "normalize": string indicating normalizaion strategy (or empty one for none)
# - "testAggDataPath": directory to put the time series list files of the test
# set which are used as inputs for forecasting evaluation (optional), the
# target time series are not contained here
# - "fullForecastAggDataPath": input + target time series of training data for
# forecasting evaluation (optional) (while "aggDataPath" only contains inputs)
# - "fullForecastTestAggDataPath": full test data (input + target) for forecasting
# evaluation (optional)
# - "testRatio" part of the dataset which should be used as test set (saved to
# separate directory, might be zero as it is optional)
savePreprocessedTimeSeries <- function(datasetName, dataPath, aggDataPath, resultsPath,
    machineNames = NULL, intervalLength = 86400, aggNamesList = NULL, aggFunctionsList = NULL,
    normalize = "", testAggDataPath = NULL, fullForecastAggDataPath = NULL,
    fullForecastTestAggDataPath = NULL, testRatio = 0) {
  cat("---Pre-processing---\n")
  startTime <- as.numeric(Sys.time())
  dir.create(resultsPath, showWarnings = FALSE, recursive = TRUE)
  if (length(list.files(aggDataPath, pattern = "\\.rds$")) == 0) {
    # Check if existing *.rds files for Spark aggregates, else pre-process CSVs
    cat("--Pre-processing CSV files from Spark (if not already done), saving as *.rds--\n")
    preprocessAggDataDir(dataPath)
    cat("--Creating time series lists for aggregations--\n")
    if (intervalLength > 0) {
      saveFixedLengthTimeSeries(datasetName = datasetName, dataPath = dataPath,
          aggDataPath = aggDataPath, resultsPath = resultsPath, machineNames = machineNames,
          intervalLength = intervalLength, aggNamesList = aggNamesList,
          aggFunctionsList = aggFunctionsList, normalize = normalize,
          testAggDataPath = testAggDataPath, fullForecastAggDataPath = fullForecastAggDataPath,
          fullForecastTestAggDataPath = fullForecastTestAggDataPath, testRatio = testRatio)
    } else {
      if (testRatio > 0) {
        cat("[WARN] Forecasting is currently not supported for time series of ",
            "variable length, corresponding parameters will be ignored.\n", sep = "")
      }
      saveVariableLengthTimeSeries(datasetName = datasetName, dataPath = dataPath,
          aggDataPath = aggDataPath, resultsPath = resultsPath, machineNames = machineNames,
          intervalLength = intervalLength, aggNamesList = aggNamesList,
          aggFunctionsList = aggFunctionsList, normalize = normalize)
    }
  } else {
    cat("--Using already pre-processed time series lists of aggregations--\n")
  }
  endTime <- as.numeric(Sys.time())
  cat("Pre-processing took ", round(endTime - startTime), "s\n", sep = "")
  return(TRUE)
}

# Entropy, internal validity, external validity table
saveStandardEvaluationTables <- function(clusterAssignmentsTable,
    dissMatrixPathList, outputPath = "./") {
  cat("---Creating cluster count table---\n")
  startTime <- as.numeric(Sys.time())
  clusterCountTable <- createClusterCountTable(clusterAssignmentsTable)
  endTime <- as.numeric(Sys.time())
  cat("Creating cluster count table took ", round(endTime - startTime), "s\n", sep = "")
  saveRDS(clusterCountTable, file = paste0(outputPath, "clusterCountTable.rds"))
  cat("---Creating entropy table---\n")
  startTime <- as.numeric(Sys.time())
  entropyTable <- createClusteringEntropyTable(clusterAssignmentsTable, normalize = TRUE)
  endTime <- as.numeric(Sys.time())
  cat("Creating entropy table took ", round(endTime - startTime), "s\n", sep = "")
  saveRDS(entropyTable, file = paste0(outputPath, "entropyTable.rds"))
  cat("---Creating internal validity table---\n")
  startTime <- as.numeric(Sys.time())
  internalCVITable <- createInternalCVITable(clusterAssignmentsTable, dissMatrixPathList)
  endTime <- as.numeric(Sys.time())
  cat("Creating internal validity table took ", round(endTime - startTime), "s\n", sep = "")
  saveRDS(internalCVITable, file = paste0(outputPath, "internalCVITable.rds"))
  cat("---Creating external validity table---\n")
  startTime <- as.numeric(Sys.time())
  externalCVITable <- createExternalCVITable(addRawResultsToClusterSummaryTable(
    clusterAssignmentsTable))
  endTime <- as.numeric(Sys.time())
  cat("Creating external validity table took ", round(endTime - startTime), "s\n", sep = "")
  saveRDS(externalCVITable, file = paste0(outputPath, "externalCVITable.rds"))
  return(TRUE)
}

# Saves multiple evaluation tables which compare how machines, sensors, days,
# days of week, weekday vs weekend are distributed over clusters
saveGroundTruthEvaluationTables <- function(clusterAssignmentsTable, tsDataBasePathList, outputPath) {
  for (i in seq_len(length(tsDataBasePathList))) {
    datasetName <- names(tsDataBasePathList)[i]
    cat("---Creating ground truth evaluation for dataset ", datasetName, "---\n", sep = "")
    allTablesStartTime <- as.numeric(Sys.time())
    tsDataDir <- tsDataBasePathList[[i]]
    cat("--Checking time series names for ground truth evaluation--\n")
    if (!haveTSListsSameNames(tsListDataDir = tsDataDir)) {
      warning("Time series lists don't have the same names, evaluation not possible.\n")
      next()
    } else {# if all time series lists have same name, we can use any of them
      tsList <- readRDS(list.files(tsDataDir, pattern = "\\.rds", full.names = TRUE)[1])
      currAssignmentsTable <- clusterAssignmentsTable[Dataset == datasetName]
      
      cat("---Creating machine ground truth table---\n")
      startTime <- as.numeric(Sys.time())
      groundTruth <- getMachineIndexVector(tsList)
      groundTruthTable <- createGroundTruthIndexTable(currAssignmentsTable, groundTruth)
      endTime <- as.numeric(Sys.time())
      cat("Creating machine ground truth table took ", round(endTime - startTime), "s\n", sep = "")
      saveRDS(groundTruthTable, file = paste0(outputPath, "groundTruthTable_", datasetName, "_machine.rds"))
      
      cat("---Creating sensor ground truth table---\n")
      startTime <- as.numeric(Sys.time())
      groundTruth <- getSensorIndexVector(tsList)
      groundTruthTable <- createGroundTruthIndexTable(currAssignmentsTable, groundTruth)
      endTime <- as.numeric(Sys.time())
      cat("Creating sensor ground truth table took ", round(endTime - startTime), "s\n", sep = "")
      saveRDS(groundTruthTable, file = paste0(outputPath, "groundTruthTable_", datasetName, "_sensor.rds"))
      
      cat("---Creating day ground truth table---\n")
      startTime <- as.numeric(Sys.time())
      groundTruth <- getDayIndexVector(tsList)
      groundTruthTable <- createGroundTruthIndexTable(currAssignmentsTable, groundTruth)
      endTime <- as.numeric(Sys.time())
      cat("Creating day ground truth table took ", round(endTime - startTime), "s\n", sep = "")
      saveRDS(groundTruthTable, file = paste0(outputPath, "groundTruthTable_", datasetName, "_day.rds"))
      
      cat("---Creating time of day ground truth table---\n")
      startTime <- as.numeric(Sys.time())
      groundTruth <- getTimeOfDayIndexVector(tsList)
      groundTruthTable <- createGroundTruthIndexTable(currAssignmentsTable, groundTruth)
      endTime <- as.numeric(Sys.time())
      cat("Creating time of day ground truth table took ", round(endTime - startTime), "s\n", sep = "")
      saveRDS(groundTruthTable, file = paste0(outputPath, "groundTruthTable_", datasetName, "_timeOfDay.rds"))
      
      cat("---Creating day of week ground truth table---\n")
      startTime <- as.numeric(Sys.time())
      groundTruth <- getDayOfWeekIndexVector(tsList)
      groundTruthTable <- createGroundTruthIndexTable(currAssignmentsTable, groundTruth)
      endTime <- as.numeric(Sys.time())
      cat("Creating day of week ground truth table took ", round(endTime - startTime), "s\n", sep = "")
      saveRDS(groundTruthTable, file = paste0(outputPath, "groundTruthTable_", datasetName, "_dayOfWeek.rds"))
      
      cat("---Creating weekday ground truth table---\n")
      startTime <- as.numeric(Sys.time())
      groundTruth <- getWeekDayVsEndIndexVector(tsList)
      groundTruthTable <- createGroundTruthIndexTable(currAssignmentsTable, groundTruth)
      endTime <- as.numeric(Sys.time())
      cat("Creating weekday ground truth table took ", round(endTime - startTime), "s\n", sep = "")
      saveRDS(groundTruthTable, file = paste0(outputPath, "groundTruthTable_", datasetName, "_weekday.rds"))
    }
    allTablesEndTime <- as.numeric(Sys.time())
    cat("Creating ground truth tables for dataset ", datasetName, " took ",
        round(allTablesEndTime - allTablesStartTime), "s\n", sep = "")
  }
  return(TRUE)
}

# Relative and absolute forecast error tables for a combination of training and
# test data, based on a clustering result and certain dissimilarities/aggregations
# (should be the same used for clustering the training data)
saveForecastEvaluationTables <- function(fullTrainDataList = NULL,
    inputTrainDataBasePathList, fullTrainDataBasePathList, trainDataPrevActList,
    trainDataAssignmentsTable, fullTestDataList = NULL, inputTestDataBasePathList,
    fullTestDataBasePathList, testDataPrevActList, testDataDissMatrixPath,
    aggList = NULL, disList = DEFAULT_DISPATCH_DISSIMILARITIES, outputPath = "./") {
  cat("---Creating forecast evaluation tables---\n")
  # Train data input + target series
  if (is.null(fullTrainDataList)) {
    fullTrainDataList <- list()
    for (datasetName in names(fullTrainDataBasePathList)) {
      fullTrainDataList[[datasetName]] <- readRDS(list.files(
        fullTrainDataBasePathList[[datasetName]], pattern = "_Raw_.*\\.rds",
        full.names = TRUE)[1])
    }
  } else {# aggregate provided dataset
      callBatchtoolsRoutine("train data (full) aggregation", FUN = aggregateAndSaveData,
                            tsListList = fullTrainDataList, aggList = aggList,
                            outputBasePathList = fullTrainDataBasePathList)
  }
  # Test data input + target series (aggregate if not already aggregated)
  if (is.null(fullTestDataList)) { # read raw test data from disk, assume aggregation
    fullTestDataList <- list()
    for (datasetName in names(fullTestDataBasePathList)) {
      fullTestDataList[[datasetName]] <- readRDS(list.files(
        fullTestDataBasePathList[[datasetName]], pattern = "_Raw_.*\\.rds",
        full.names = TRUE)[1])
    }
    inputTestDataList <- lapply(names(fullTestDataList), function(datasetName) {
      fullTestDataList[[datasetName]][testDataPrevActList[[datasetName]]$prev]
    })
    names(inputTestDataList) <- names(fullTestDataList)
  } else {# aggregate the provided dataset
    callBatchtoolsRoutine("test data (full) aggregation", FUN = aggregateAndSaveData,
        tsListList = fullTestDataList, aggList = aggList,
        outputBasePathList = fullTestDataBasePathList)
    inputTestDataList <- lapply(names(fullTestDataList), function(datasetName) {
      fullTestDataList[[datasetName]][testDataPrevActList[[datasetName]]$prev]
    })
    names(inputTestDataList) <- names(fullTestDataList)
    callBatchtoolsRoutine("test data (input) aggregation", FUN = aggregateAndSaveData,
        tsListList = inputTestDataList, aggList = aggList,
        outputBasePathList = inputTestDataBasePathList)
  }
  # Assign test data to clusters
  clusterAssignmentsTable <- formatAndSortClusterSummaryTable(trainDataAssignmentsTable)
  testDataAssignmentResult <- callBatchtoolsRoutine("test data assignment",
      FUN = assignDatasetsToClusters.dissParallelized,
      newDataBasePathList = inputTestDataBasePathList,
      oldDataBasePathList = inputTrainDataBasePathList,
      clusterAssignmentsTable = clusterAssignmentsTable, disList = disList,
      numNeighbors = 1, dissMatrixPath = testDataDissMatrixPath)
  testDataAssignmentTimeTable <- extractExecutionTimeTable(testDataAssignmentResult)
  saveRDS(testDataAssignmentResult, file = paste0(outputPath, "testDataAssignmentResult.rds"))
  saveRDS(testDataAssignmentTimeTable, file = paste0(outputPath, "testDataAssignmentTimeTable.rds"))
  # Forecast - raw data (combined with assignments from aggregation levels)
  cat("--Creating forecast evaluation table on base level--\n")
  startTime <- as.numeric(Sys.time())
  baseForecastTable <- createBaseForecastErrorTable(
    newDatasetList = fullTestDataList, oldDatasetList = fullTrainDataList,
    newDataPrevActList = testDataPrevActList, oldDataPrevActList = trainDataPrevActList,
    newDataAssignmentsTable = testDataAssignmentResult,
    oldDataAssignmentsTable = trainDataAssignmentsTable,
    dissMatrixPaths = testDataDissMatrixPath, dissWeightConversion = "gaussian"
  )
  endTime <- as.numeric(Sys.time())
  cat("Creating base forecast evaluation table took ", round(endTime - startTime), "s\n", sep = "")
  saveRDS(baseForecastTable, file = paste0(outputPath, "baseForecastTable.rds"))
  # Forecast - aggregated data (data used for the clustering on each level)
  cat("--Creating relative forecast evaluation table--\n")
  startTime <- as.numeric(Sys.time())
  relativeForecastTable <- createRelativeForecastErrorTable(
    newDatasetBasePathList = fullTestDataBasePathList,
    oldDatasetBasePathList = fullTrainDataBasePathList,
    newDataPrevActList = testDataPrevActList, oldDataPrevActList = trainDataPrevActList,
    newDataAssignmentsTable = testDataAssignmentResult,
    oldDataAssignmentsTable = trainDataAssignmentsTable,
    dissMatrixPaths = testDataDissMatrixPath, dissWeightConversion = "gaussian",
    forecastAggregation = "mean"
  )
  endTime <- as.numeric(Sys.time())
  cat("Creating relative forecast evaluation table took ", round(endTime - startTime), "s\n", sep = "")
  saveRDS(relativeForecastTable, file = paste0(outputPath, "relativeForecastTable.rds"))
  return(TRUE)
}

# Runs a full impact analysis (without forecasting evaluation) with the steps
# 1) Aggregation (optional, only if dataset directly passed)
# 2) Clustering (plots are optional)
# 3) Entropy / internal / external evaluation
# 4) Ground truth evaluation (properties in clusters) (optional)
# Results are saved as *.rds files
#
# Parameters:
# - "datasetList": list of list of time series to be clustered; might be NULL if
# aggregated datasets are already located under "aggDataBasePathList" and forecasting
# evaluation won't be needed ("testDataBasePathList" also NULL)
# - "aggDataBasePathList": directory where the aggregations of the datasets are saved
# - "dissMatrixPathList": directory where dissimilarity matrices are saved
# - "tsPlotPathList": list of directories where the datasets are plotted (optional);
# if provided, the plots will be copied according to the clustering results
# - "clusPlotBasePathList": list of directories for cluster plots (members, centroid)
# and target of time series plot copies (previous parameter) (optional)
# - "resultsPath": output directory for *.rds file created by analysis
# - "aggList": list of aggregation methods (a default is provided, parameter
# is only used if data are not pre-aggregated or forecasts should be made)
# - "disList": list of dissimilarity methods (a default is provided)
# - "clusList": list of clustering methods (a default is provided)
# - "groundTruth": run ground truth evaluation (requires special naming scheme of
# time series)
runImpactAnalysis <- function(datasetList = NULL, aggDataBasePathList,
    dissMatrixPathList, tsPlotPathList = list(), clusPlotBasePathList = list(),
    resultsPath = "./", aggList = DEFAULT_AGGREGATIONS,
    disList = DEFAULT_DISPATCH_DISSIMILARITIES, clusList = DEFAULT_CLUSTERINGS,
    groundTruth = TRUE) {
  dir.create(resultsPath, showWarnings = FALSE, recursive = TRUE)
  # Aggregation
  if (!is.null(datasetList)) {
    callBatchtoolsRoutine("data aggregation", FUN = aggregateAndSaveData,
        tsListList = datasetList, aggList = aggList, outputBasePathList = aggDataBasePathList)
  }
  # Clustering
  clusteringResult <- callBatchtoolsRoutine("clustering",
      FUN = analyzeAggregationImpact.preagg, aggBasePathList = aggDataBasePathList,
      disList = disList, clusList = clusList, dissMatrixPathList = dissMatrixPathList,
      tsPlotPathList = tsPlotPathList, clusPlotBasePathList = clusPlotBasePathList)
  clusteringTimeTable <- extractExecutionTimeTable(clusteringResult)
  saveRDS(clusteringResult, paste0(resultsPath, "clusteringResult.rds"))
  saveRDS(clusteringTimeTable, paste0(resultsPath, "clusteringTimeTable.rds"))
  # Evaluation - Entropy, Internal, External
  saveStandardEvaluationTables(clusteringResult, dissMatrixPathList, resultsPath)
  # Evaluation - ground truth (properties describing time series in clusters)
  if (groundTruth) {
    saveGroundTruthEvaluationTables(clusteringResult, aggDataBasePathList, resultsPath)
  }
  return(TRUE)
}

# Parses intput parameters and saves them in a file
parseInputParams <- function(outputFile = "params.txt") {
  params <- parse_args(object = OptionParser(option_list = list(
    make_option("--dataPath", action = "store", type = "character", default = "./"),
    make_option("--plotPath", action = "store", type = "character", default = ""),
    make_option("--datasetName", action = "store", type = "character", default = "MyDataset"),
    make_option("--normalize", action = "store", type = "character", default = ""),
    make_option("--testRatio", action = "store", type = "numeric", default = 0),
    make_option("--aggregations", action = "store", type = "character",
                default = "mean|median|min|max|stddev|skewness|kurtosis"),
    make_option("--machines", action = "store", type = "character",
                default = ""),
    make_option("--intervalLength", action = "store", type = "numeric",
                default = 86400)
  ), add_help_option = FALSE))
  params <- checkParameters(params)
  outputDir <- paste0(params$dataPath, "results/")
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  sink(paste0(outputDir, outputFile))
  cat("Time: ", format(Sys.time()), "\n", sep = "")
  cat("Commit: ", system("git log --oneline -1", intern = TRUE), "\n", sep = "")
  cat("Parameters:\n")
  cat(paste0("(", seq_len(length(params)), ") ", names(params), ": ", params, "\n"), sep = "")
  sink()
  saveRDS(params, file = paste0(outputDir, "params.rds"))
  return(params)
}

# Checks the values of input parameters to the batch script and either stops if
# values are invalid or creates a warning and fixes them (the changed parameter
# list is returned).
checkParameters <- function(params) {
  if (substr(params$dataPath, start = nchar(params$dataPath),
             stop = nchar(params$dataPath)) != "/") {
    params$dataPath <- paste0(params$dataPath, "/")
    cat("[INFO] --dataPath did not end with a slash; we added one.")
  }
  if (!dir.exists(params$dataPath)) {
    stop("--dataPath is not an existing directory.")
  }
  if (params$plotPath != "" && substr(params$plotPath, start = nchar(params$plotPath),
                                      stop = nchar(params$plotPath)) != "/") {
    params$plotPath <- paste0(params$plotPath, "/")
    cat("[INFO] --plotPath did not end with a slash; we added one.")
  } # no existence check, can also be created during analysis
  if (grepl(pattern = "_", x = params$datasetName)) {
    stop("Because of internal naming conventions, --datasetName must not contain
         an undercore.")
  }
  if (!(params$normalize %in% c("", "zscore", "minmax"))) {
    stop("--normalize has a value which is not supported.")
  }
  if (params$testRatio >= 1 || params$testRatio < 0) {
    stop("--testRatio has to be a number in [0,1).")
  }
  if (params$intervalLength <= 0) {
    cat("[INFO] --intervalLength <= 0 will create time series of variable length.\n")
    if (params$intervalLength == 0 & !file.exists(paste0(params$dataPath,
        "results/timeSeriesOverviewTableInput.rds"))) {
      stop(paste0("Using --intervalLength of 0 requires a file 'timeSeriesOverviewTableInput.rds' ",
                    "placed at a directory 'results' under the --dataPath."))
    }
    if (params$testRatio > 0) {
      params$testRatio <- 0
      cat("[INFO] Forecast evaluation is disabled for variable length time series.\n")
    }
  } else if (86400 / params$intervalLength != floor(86400 / params$intervalLength)) {
    stop("--intervalLength must divide the seconds of a day (86400) without remainder.")
  }
  return(params)
}
