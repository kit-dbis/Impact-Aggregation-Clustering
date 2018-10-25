library(data.table)
library(FastTSDistances)
library(zoo)

### Functions ###

# Checks all *.rds files in a directory which represent lists if they all have
# the same names.
haveTSListsSameNames <- function(tsListDataDir = "./") {
  tsNames <- character()
  firstTSListFile <- ""
  for (fileName in list.files(tsListDataDir, pattern = "\\.rds$", full.names = TRUE)) {
    tsList <- readRDS(fileName)
    if (is.list(tsList)) {
      if (length(tsNames) == 0) {
        tsNames <- names(tsList)
        firstTSListFile <- fileName
      } else if (length(tsList) != length(tsNames)) {
        warning(paste0("File ", fileName, " has a different number of time series",
                       " than ", firstTSListFile, "."))
        return(FALSE)
      } else if (any(names(tsList) != tsNames)) {
        warning(paste0("File ", fileName, " has different names for its lists than ",
                        firstTSListFile, "."))
        return(FALSE)
      }
    }
  }
  if (length(tsNames) == 0) {
    warning("No files with time series lists found.")
  }
  return(TRUE)
}

# Handles missing values in a data.table per machine and sensor for a certain column;
# modification of data.table in place (!!)
# "strategy": "linear" (default), "cf" for carry-forward, "cb" for carry-backward;
# leading/following NAs are replaced by cb/cf (so series with at least one non-NA value
# are definitely free of NAs after this procedure), pure NA series remain as-is
handleMissingValues <- function(sensorDataTable, strategy = "", colName = "value") {
  setorder(sensorDataTable, ts)
  if (strategy == "cf") {
    sensorDataTable[, (colName) := na.locf(na.locf(get(colName), na.rm = FALSE),
        fromLast = TRUE, na.rm = FALSE), by = .(machine_name, sensor_name)]
  } else if (strategy == "cb") {
    sensorDataTable[, (colName) := na.locf(na.locf(get(colName), na.rm = FALSE,
        fromLast = TRUE), na.rm = FALSE), by = .(machine_name, sensor_name)]
  } else {
    sensorDataTable[, (colName) := approxImproved(x = ts, y = get(colName), xout = ts)$y,
          by = .(machine_name, sensor_name)]
    sensorDataTable[, (colName) := na.locf(na.locf(get(colName), na.rm = FALSE),
        fromLast = TRUE, na.rm = FALSE), by = .(machine_name, sensor_name)] # NAs at beginning and end
  }
  return(sensorDataTable)
}

# Can handle constant and NA vectors
approxImproved <- function(x, y, xout) {
  if (uniqueN(y) == 1) {
    return(list(x = x, y = y))
  } else {
    return(approx(x, y, xout))
  }
}

# Removes all time series which have less than one distinct non-NA value from a
# data.table or a list of time series
removeConstantSeries <- function(sensorData) {
  if (is.data.table(sensorData)) {
    uniqueNTable <- sensorData[, .(DistinctCount = uniqueN(value, na.rm = TRUE)),
                               by = .(machine_name, sensor_name)]
    return(merge(sensorData, uniqueNTable[DistinctCount > 1]))
  } else {
    return(sensorData[sapply(sensorData, function(ts) uniqueN(ts, na.rm = TRUE) > 1)])
  }
}

# Returns the statistical mode parameter (only the first one if multi-modal)
# copied from https://stackoverflow.com/a/25635740
modeOfData <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Segments a numerical vector by ignoring parts which equal a certain value
# "splitValue": value which will be used for segmentation; mode of series if
# not provided
# "gapTolerance": subsequent occurences of the split value up to this number will
# not cause finishing the current segment
# "numTolerance": numerical threshold for equality
# returns: a list of numeric vectors representing segments of the original series
segmentByValue <- function(series, splitValue = modeOfData(series), gapTolerance = 0,
                           numTolerance = 0) {
  isEqualMode <- abs(series - splitValue) <= numTolerance
  # Determine start and end with help of mode
  segmentStartIdx <- which(!isEqualMode & shift(isEqualMode, type = "lag"))
  if (!isEqualMode[1]) {
    segmentStartIdx <- c(1, segmentStartIdx)
  }
  segmentEndIdx <- which(!isEqualMode & shift(isEqualMode, type = "lead"))
  if (!isEqualMode[length(isEqualMode)]) {
    segmentEndIdx <- c(segmentEndIdx, length(series))
  }
  # Analyze if there are very short sequences equaling the mode and ignore them
  betweenSegmentLength <- (segmentStartIdx[-1] - segmentEndIdx[-length(segmentEndIdx)]) - 1
  segmentStartIdx <- segmentStartIdx[c(TRUE, betweenSegmentLength > gapTolerance)]
  segmentEndIdx <- segmentEndIdx[c(betweenSegmentLength > gapTolerance, TRUE)]
  # For each segment start: extract points belonging to segment
  return(lapply(1:length(segmentStartIdx), function(segmentIdx)
    series[segmentStartIdx[segmentIdx]:segmentEndIdx[segmentIdx]]))
}

# Takes 
# - a numerical vector and converts it to a list of segments with the
# "targetLength" (apart from the last entry which can be shorter)
# - a data.table and converts its "value" and "ts" attribute to a list of segments
# based on the "ts" attribute and "targetLength" (which will be interpreted as interval
# length in seconds); the timestamps will be used as names
segmentToLength <- function(series, targetLength = 100) {
  if (is.numeric(series)) {
    valuesList <- lapply(1:(ceiling(length(series) / targetLength)), function(i) {
      return(series[((i - 1)*targetLength + 1):min(i*targetLength, length(series))])
    })
  } else if (is.data.table(series)) {
    minTimeStamp <- series[, as.numeric(min(ts))]
    splitLevels <- series[, floor((as.numeric(ts) - minTimeStamp) / targetLength)]
    valuesList <- split(series$value, splitLevels)
    tsList <- split(series$ts, splitLevels)
    for (i in 1:length(valuesList)) {
      names(valuesList[[i]]) <- tsList[[i]]
    }
  } else {
    return(NA)
  }
  names(valuesList) <- sapply(1:length(valuesList), function(i) paste0("Segment_", i))
  return(valuesList)
}

# Conversion from factor/string to list(POSIXct, numeric, character)
# [timestamp, numeric timestamp, year-month-day]
# maxSimultaneousTimestampConversion due to inefficient memory usage in as.POSICxt()
# -> data need to be split and processed in splits, results merged
convertTimestamps.aggData <- function(oldTimestamps, maxSimultaneousTimestampConversion = 1e6) {
  oldTimestamps <- as.character(oldTimestamps)
  newPOSIX <- as.POSIXct(numeric(length(oldTimestamps)), origin = "1970-01-01", tz = "GMT")
  newNumeric <- numeric(length = length(oldTimestamps))
  newDayChar <- character(length = length(oldTimestamps))
  newSecOfDay <- numeric(length = length(oldTimestamps))
  dataSplitCount <- ceiling(length(oldTimestamps) / maxSimultaneousTimestampConversion)
  gc()
  for (i in 1:dataSplitCount) {
    print(paste0("Converting timestamps - split ", i, "/", dataSplitCount, " ..."))
    startIndex <- 1 + (i - 1) * maxSimultaneousTimestampConversion
    endIndex <- min(length(oldTimestamps), i * maxSimultaneousTimestampConversion)
    # old timestamp like "2017-02-01T23:00:00.000+01:00", we ignore time zone
    newPOSIX[startIndex:endIndex] <- as.POSIXct(oldTimestamps[startIndex:endIndex],
        format = "%Y-%m-%dT%H:%M:%S", tz = "GMT") # use UTC which has no daylight-saving time
    newNumeric[startIndex:endIndex] <- as.numeric(newPOSIX[startIndex:endIndex])
    newDayChar[startIndex:endIndex] <- substr(oldTimestamps[startIndex:endIndex], 1, 10)
    # Seconds of day in local time (ignoring leap seconds, summer time)
    newSecOfDay[startIndex:endIndex] <-
      as.integer(substr(oldTimestamps[startIndex:endIndex], 12, 13)) * 3600 +
      as.integer(substr(oldTimestamps[startIndex:endIndex], 15, 16)) * 60 +
      as.integer(substr(oldTimestamps[startIndex:endIndex], 18, 19))
    gc()
  }
  return(list(newPOSIX, newNumeric, newDayChar, newSecOfDay))
}

# Creates "numberOfTS" equi-spaced timestamps, starting with "from" and considering the
# "samplingInterval" in seconds. Returns a list with the timestamps as POSIXct,
# timestamps as numeric, the day as character and the second of day as numeric
createTimestamps.aggData <- function(from, samplingInterval, numberOfTS,
                                     maxSimultaneousTimestampConversion = 1e6) {
  ts <- seq(from = from, by = samplingInterval, length.out = numberOfTS)
  ts_numeric <- numeric(length = length(ts))
  day <- character(length = length(ts))
  secOfDay <- numeric(length = length(ts))
  dataSplitCount <- ceiling(length(ts) / maxSimultaneousTimestampConversion)
  gc()
  for (i in 1:dataSplitCount) {
    cat("Converting timestamps - split ", i, "/", dataSplitCount, "...\n", sep = "")
    startIndex <- 1 + (i - 1) * maxSimultaneousTimestampConversion
    endIndex <- min(length(ts), i * maxSimultaneousTimestampConversion)
    ts_numeric[startIndex:endIndex] <- as.numeric(ts[startIndex:endIndex])
    day[startIndex:endIndex] <- format(ts[startIndex:endIndex], format = "%Y-%m-%d")
    secOfDay[startIndex:endIndex] <- ts_numeric[startIndex:endIndex] -
      as.numeric(trunc(ts[startIndex:endIndex], units = "days"))
    gc()
  }
  return(list(ts, ts_numeric, day, secOfDay))
}

# Read, make data.table, rename attributes, convert timestamps, sort, save as RDS
loadAndPreprocessAggData <- function(inputFilePath, outputFilePath = NULL) {
  print("Reading input ...")
  aggData <- data.table(read.csv(file = inputFilePath))
  print("Pre-processing ...")
  aggData <- aggData[, -"date_end"]
  setnames(aggData, old = "date_start", new = "ts")
  # Simultaneous conversion of all timestamps would consume too much memory for
  # large datasets - therefore some more code here
  aggData[, c("ts", "ts_numeric", "day", "secOfDay") := convertTimestamps.aggData(ts)]
  gc()
  # Move new time columns from end to position directly after old time column
  # (which is at position 3)
  setcolorder(aggData, c(1:3, (ncol(aggData) - 2):ncol(aggData), 4:(ncol(aggData) - 3)))
  print("Sorting")
  setorder(aggData, ts_numeric)
  if (!is.null(outputFilePath)) {
    print("Writing output ...")
    saveRDS(aggData, file = outputFilePath)
  }
  print("Done.")
  return(aggData)
}

# Pre-processes all aggregated datasets from CSV files in a directory and writes
# the results into CSVs with the same file name
preprocessAggDataDir <- function(dirPath = "./") {
  for (fileName in list.files(path = dirPath, pattern = "\\.csv$")) {
    cat("##### File: ", fileName, "######\n")
    outputFilePath <- paste0(dirPath, gsub("csv$", "rds", fileName))
    if (file.exists(outputFilePath)) {
      cat("Already pre-processed.\n")
    } else {
      loadAndPreprocessAggData(inputFilePath = paste0(dirPath, fileName),
                               outputFilePath = outputFilePath)
    }
  }
  return(TRUE)
}

# Reads in all RDS files in a dir (which should be the result of "preprocessAggDataDir()")
# and limits the datasets to certain machine, sensor, day combinations which are
# passed as parameter. Writes the resulting files to the same dir with a suffix.
reduceAggDataDirToValidDays <- function(dirPath, machineSensorDayTable) {
  inputRDSFiles <- list.files(path = dirPath, pattern = "\\.rds$", full.names = TRUE)
  inputRDSFiles <- inputRDSFiles[!endsWith(inputRDSFiles, "_valid.rds")]
  cat("-----Limiting datasets to valid days.-----\n")
  for (inputFilePath in inputRDSFiles) {
    cat("File:", inputFilePath, "\n")
    outputFilePath <- gsub("\\.rds$", "_valid.rds", inputFilePath)
    if (file.exists(outputFilePath)) {
      cat("Already pre-processed.\n")
    } else {
      dataset <- readRDS(inputFilePath)
      dataset <- merge(dataset, machineSensorDayTable)
      setorder(dataset, ts_numeric)
      saveRDS(dataset, outputFilePath)
    }
  }
}

# Runs the standard CSV -> RDS pre-processing and reduces the dataset afterwards
# by limiting it to certain days; "validDay" function should return all days
# which are valid (called per machine-sensor)
preprocessAndReduceAggDataDir <- function(dataPath, validDayFunction) {
  # Standard CSV -> RDS conversion
  preprocessAggDataDir(dataPath)
  # Find valid combinations with biggest file
  inputRDSFiles <- list.files(path = dataPath, pattern = "\\.rds$", full.names = TRUE)
  inputRDSFiles <- inputRDSFiles[!endsWith(inputRDSFiles, "_valid.rds")]
  datasetFile <- inputRDSFiles[which.max(file.size(inputRDSFiles))]
  dataset <- readRDS(datasetFile)
  validMachineDayCombis <- dataset[, .(day = do.call(validDayFunction, list(.SD))),
                                   by = .(machine_name, sensor_name)]
  rm(dataset)
  # Reduce
  reduceAggDataDirToValidDays(dirPath = dataPath,
                              machineSensorDayTable = validMachineDayCombis)
  return(TRUE)
}

# Prints some diagnosis statistics for a aggregates dataset created by Spark
# (NAs, missing values); should only be used after calling "loadAndPreprocess"
diagnoseAggData <- function(aggDataTable) {
  print("---NAs per column---")
  print(sapply(aggDataTable, function(column) sum(is.na(column))))
  print("---Values per machine/sensor---")
  print(aggDataTable[, .N, keyby = .(machine_name, sensor_name)])
  print("---Values per day---")
  valuesPerDay <- aggDataTable[, .(ValuesCount = .N), keyby = day]
  daysWithValuesCountCount <- valuesPerDay[, .N, keyby = ValuesCount]
  print(daysWithValuesCountCount)
  print("---Days with exceptional number of values---")
  print(valuesPerDay[ValuesCount != daysWithValuesCountCount[N == max(N), ValuesCount]])
  print("---Machines/sensors and days without any values---")
  # expected: all combinations of machines/sensors and days should exist in dataset
  machineSensorDayTable <- merge(x = aggDataTable[, 1, by = .(machine_name, sensor_name)],
                                 y = data.table(V1 = 1, day = aggDataTable[, unique(day)]),
                                 by = "V1", allow.cartesian = TRUE)[, -"V1"]
  machineSensorDayTable <- merge(x = machineSensorDayTable,
      y = aggDataTable[, .N, keyby = .(machine_name, sensor_name, day)],
      by = colnames(machineSensorDayTable), all.x = TRUE)
  print(machineSensorDayTable[is.na(N), -"N"])
  print("--Machine/sensors and days with less values than expected---")
  expectedCount <- machineSensorDayTable[, modeOfData(N)]
  print(machineSensorDayTable[!is.na(N) & N != expectedCount][order(day)])
  return(TRUE)
}

# removes all combinations of machines and sensors and days which do not have
# at least "expectedCount" values (absolute or a fraction of the most frequent
# number of values of any machine-day combination)
removeSparseMachineDaySensorCombis <- function(aggDataTable, expectedCount) {
  machineSensorDayCounts <- aggDataTable[, .N, by = .(machine_name, sensor_name, day)]
  if (expectedCount <= 1) {
    expectedCount <- expectedCount * machineSensorDayCounts[, modeOfData(N)]
  }
  result <- merge(x = machineSensorDayCounts[N >= expectedCount][, -"N"],
                  y = aggDataTable, all.x = TRUE)
  setcolorder(result, c(1, 3:5, 2, 6:ncol(aggDataTable))) # move "day" column
  return(result)
}

# Make sure sensor values appear in certain interval
#
# - "samplingInterval": value which can be used for creating a POSIXct sequence; if not
# provided, most frequent sampling interval in dataset is used
# - "from"/"to": either as POSIXct or min/max in dataset are taken if NA
# - "strategy": "linear" (default), "cf" for carry-forward, "cb" for carry-backward;
# leading/following NAs are replaced by cb/cf (so series with at least one non-NA value
# are definitely free of NAs after this procedure), pure NA series remain as-is
# - "onlyExistingCombinations" decides if only machine-sensor-day-timeOfDay combinations
# from the data should be interpolated or really all possible ones (intervals with no
# measurements will be completely interpolated in that case);
# - "intervalLength": if positive: number of seconds represented by one time
# series, default is one day; if negative: variable length time series should
# be possible, therefore no fixed time intervals and "onlyEistingCombinations"
# is ignored
interpolateAllAggData <- function(aggDataTable, from = NA, to = NA, samplingInterval = NA,
    strategy = "linear", onlyExistingCombinations = TRUE, intervalLength = 86400) {
  if (is.na(from)) {
    from <- aggDataTable[, min(ts)]
  }
  if (is.na(to)) {
    to <- aggDataTable[, max(ts)]
  }
  if (is.na(samplingInterval)) { # take most frequent difference between timestamps in dataset
    samplingInterval <- modeOfData(aggDataTable[order(ts)][, diff(ts_numeric),
        by = .(machine_name, sensor_name)]$V1)
  }
  timestampsTable <- data.table(ts = seq(from = from, to = to, by = samplingInterval))
  timestampsTable[, ts_numeric := as.numeric(ts)]
  if (intervalLength > 0) {
    # These attributes only matter if fixed-interval time series are created
    timestampsTable[, day := format(ts, format = "%Y-%m-%d")]
    timestampsTable[, secOfDay := as.integer(format(ts, format = "%H")) * 3600 +
                      as.integer(format(ts, format = "%M")) * 60 +
                      as.integer(format(ts, format = "%S"))]
    timestampsTable <- addTimeOfDay(timestampsTable, intervalLength = intervalLength)
  }
  gc() # because of POSIXct operation above
  if (onlyExistingCombinations & intervalLength > 0) {
    machinesAndDays <- aggDataTable[, 1, by = .(machine_name, sensor_name,
                                                day, timeOfDay)][, -"V1"]
    result <- merge(machinesAndDays, timestampsTable, by = c("day", "timeOfDay"),
                    allow.cartesian = TRUE)
    result <- merge(x = result, y = aggDataTable, all.x = TRUE, all.y = TRUE,
                    by = c("machine_name", "sensor_name", "ts", "ts_numeric", "day",
                           "secOfDay", "timeOfDay"))
  } else {# all timestamps for all machine-sensor combinations
    timestampsTable[, V1 := 1] # merge dummy
    result <- merge(aggDataTable[, 1, by = .(machine_name, sensor_name)],
                    timestampsTable, allow.cartesian = TRUE)[, -"V1"]
    result <- merge(x = result, y = aggDataTable, all.x = TRUE, all.y = TRUE)
  }
  setorder(result, ts)
  # For all aggregate columns -> interpolate/replace
  for (aggColumn in setdiff(colnames(result), c("machine_name", "sensor_name", "ts",
        "ts_numeric", "day", "secOfDay", "timeOfDay"))) {
    result[, (aggColumn) := as.numeric(get(aggColumn))]
    handleMissingValues(result, strategy = strategy, colName = aggColumn)
  }
  result <- result[ts_numeric %in% timestampsTable$ts_numeric] # use only new timestamps
  return(result)
}

# In-place ordering of a table containing sensor data
orderForAggDataTimeSeriesListConversion <- function(aggDataTable) {
  if ("ts_numeric" %in% colnames(aggDataTable)) {
    setorder(aggDataTable, machine_name, sensor_name, ts_numeric)
  } else if ("ts" %in% colnames(aggDataTable)) {
    setorder(aggDataTable, machine_name, sensor_name, ts)
  } else {
    setorder(aggDataTable, machine_name, sensor_name, day, timeOfDay)
  }
  return(TRUE)
}

# Converts a data.table to a list of time series based on the "splitLevels",
# using one or multiple attributes as values of the time series and the timestamps
# as names
convertToTimeSeriesList.aggData <- function(aggDataTable, splitLevels, attributeNames,
                                            normalize = "") {
  if (length(attributeNames) == 1) {# univariate time series
    valuesList <- split(x = aggDataTable[, get(attributeNames)], f = splitLevels)
  } else {# multivariate time series
    valuesList <- lapply(split(x = aggDataTable[, mget(attributeNames)],
                               f = splitLevels), as.matrix)
  }
  timestampsList <- split(x = aggDataTable[, ts], f = splitLevels)
  for (i in 1:length(valuesList)) {
    names(valuesList[[i]]) <- timestampsList[[i]]
  }
  if (normalize == "zscore") {
    valuesList <- lapply(valuesList, znormalize)
  } else if (normalize == "minmax") {
    valuesList <- lapply(valuesList, minMaxNormalize)
  }
  return(valuesList)
}

# Takes data.table of sensor values and creates a list with one time series per
# machine and sensor and day and timeOfDay (interval);
# "attributeNames" denote the columns forming the time  series attributes, might
# be just one (univariate time series, numeric vectors will be created) or
# multiple ones (mutivariate time series, numeric matrices will be created)
convertToFixedIntervalTimeSeriesList <- function(aggDataTable, attributeNames = "mean",
                                                 normalize = "") {
  if (!"timeOfDay" %in% names(aggDataTable)) { # make full day time series
    aggDataTable <- copy(aggDataTable)
    aggDataTable[, timeOfDay := 1]
  }
  aggDataTable <- aggDataTable[, mget(c("machine_name", "sensor_name", "day",
                                        "timeOfDay", "ts", attributeNames))]
  # Ordering important if this routines is called for multiple aggregated datasets
  # which might have different input orderings (should not happen, but make sure anyway)
  orderForAggDataTimeSeriesListConversion(aggDataTable)
  splitLevels <- aggDataTable[, paste0(machine_name, "-", sensor_name, "-",
                                       gsub("-", "_", day), "-", timeOfDay)]
  return(convertToTimeSeriesList.aggData(aggDataTable = aggDataTable,
      splitLevels = splitLevels, attributeNames = attributeNames, normalize = normalize))
}

# Takes data.table of sensor values and creates a list with one time series per
# row in "timeSeriesCombinations" (which should contain machine_name, sensor_name,
# start time and end time of the series to be created);
# - "attributeNames" denote the columns forming the time  series attributes, might
# be just one (univariate time series, numeric vectors will be created) or
# multiple ones (mutivariate time series, numeric matrices will be created)
# - "normalize" enables z-scoring or another normalization per time series
convertToVariableLengthTimeSeriesList <- function(aggDataTable, timeSeriesCombinations,
    attributeNames = "mean", normalize = "") {
  aggDataTable <- aggDataTable[, mget(c("machine_name", "sensor_name", "ts", "ts_numeric",
                                        attributeNames))]
  # Ordering important if this routines is called for multiple aggregated datasets
  # which might have different input orderings (should not happen, but make sure anyway)
  orderForAggDataTimeSeriesListConversion(aggDataTable)
  samplingInterval <- modeOfData(aggDataTable[, diff(ts_numeric)])
  gc() # mode operation might consume some memory (and not free it, whyever)
  setnames(aggDataTable, "ts_numeric", "ts_numeric_start")
  aggDataTable[, ts_numeric_end := ts_numeric_start + samplingInterval]
  # Copy to prevent modyfing the original object
  timeSeriesCombinations <- copy(timeSeriesCombinations)
  # Reduce intervals by 1 sec at begin and end because foverlaps will count
  # matches if the aggregation intervals just ends at the time when the time
  # series interval starts (which is only one point overlap)
  timeSeriesCombinations[, startTime_numeric := startTime_numeric + 1]
  timeSeriesCombinations[, endTime_numeric := endTime_numeric - 1]
  setkey(timeSeriesCombinations, machine_name, sensor_name, startTime_numeric,
         endTime_numeric) # necessary for foverlaps
  # Merge/intersection with "foverlaps" (retain only rows in aggDataTable which
  # overlap with at least one row in timeSeriesCombinations)
  aggDataTable <- foverlaps(aggDataTable, timeSeriesCombinations, type = "any", nomatch = 0,
      by.x = c("machine_name", "sensor_name", "ts_numeric_start", "ts_numeric_end"))
  splitLevels <- aggDataTable[, paste0(machine_name, "-", sensor_name, "-",
                                       gsub("-", "_", startDay), "-", startTimeOfDay)]
  return(convertToTimeSeriesList.aggData(aggDataTable = aggDataTable,
      splitLevels = splitLevels, attributeNames = attributeNames, normalize = normalize))
}

# Adds a new column denoting the current time interval of the day based on the
# intervalLength (default: one interval per day)
addTimeOfDay <- function(aggDataTable, intervalLength = 86400) {
  return(aggDataTable[, timeOfDay := as.integer(secOfDay / intervalLength + 1)])
}

# Reads all RDS files from the "inputDataDir", assuming that they are pre-processed
# tables ("preprocessAggDataDir()") containing aggregated data. For each dataset,
# the relative frequency of values per machine and sensor and interval of day
# (default: full day intervals) is counted and this is merged over the datasets.
# An additional column is added which contains the minimum frequency of values over
# all datasets. This might be used to decide which combinations of machines and
# sensors and days/intervals of days should be chosen for creating time series.
summarizeValueCountForAggDataDir <- function(inputDataDir = "./", intervalLength = 86400) {
  result <- NULL
  for (fileName in list.files(path = inputDataDir, pattern = "rds$")) {
    aggDataTable <- readRDS(file = paste0(inputDataDir, fileName))
    aggDataTable <- addTimeOfDay(aggDataTable, intervalLength = intervalLength)
    machineDaySensorCounts <- aggDataTable[, .N,
          by = .(machine_name, sensor_name, day, timeOfDay)]
    maxCount <- machineDaySensorCounts[, max(N)]
    if (maxCount <= 1) {
      cat("[WARN] The dataset file \"" , fileName, "\" has at most 1 point per ",
          "interval. This might result in problems with some dissimilarities and ",
          "therefore be an indicator to remove this aggregation level.\n", sep = "")
    }
    modeCount <- machineDaySensorCounts[, modeOfData(N)]
    if (maxCount != modeCount) {
      cat("[WARN]  The maximum number of values per time interval (", maxCount,
          ") is not equal to the most frequent number of values per time interval ",
          "(", modeCount, ") for the file \"", fileName, "\". We take the mode as ",
          "expected number of values.\n", sep = "")
    }
    machineDaySensorCounts[, N := N / modeCount]
    setnames(machineDaySensorCounts, old = "N", new = paste0("NRatio.",
        regmatches(fileName, regexpr("[0-9]+_[a-z]+", fileName))))
    if (is.null(result)) {
      result <- machineDaySensorCounts
    } else {
      result <- merge(result, machineDaySensorCounts, all.x = TRUE, all.y = TRUE)
    }
  }
  result[is.na(result)] <- 0
  result[, NRatio := do.call(pmin, .SD),
         .SDcols = grep("NRatio", colnames(result), value = TRUE)]
  return(result)
}

# Reads all RDS files from the "inputDataDir", assuming that they are pre-processed
# tables ("preprocessAggDataDir()") containing aggregated data. For each dataset,
# it is evaluated if there is only one distinct value per machine and sensor and
# interval of the day (default: there is only one interval of 86400s per day, but
# you might also pass smaller numbers [which should be a divisor of 86400]).
# The result might be used to exclude certain combinations of machines and days/
# intervals of days when creating time series.
summarizeConstantSeriesForAggDataDir <- function(inputDataDir = "./", attribute,
                                                 intervalLength = 86400) {
  result <- NULL
  for (fileName in list.files(path = inputDataDir, pattern = "rds$")) {
    aggDataTable <- readRDS(file = paste0(inputDataDir, fileName))
    aggDataTable <- addTimeOfDay(aggDataTable, intervalLength = intervalLength)
    machineDaySensorCounts <- aggDataTable[, uniqueN(get(attribute)) == 1,
        by = .(machine_name, sensor_name, day, timeOfDay)]
    setnames(machineDaySensorCounts, old = "V1", new = paste0("constant.",
        regmatches(fileName, regexpr("[0-9]+_[a-z]+", fileName))))
    if (is.null(result)) {
      result <- machineDaySensorCounts
    } else {
      result <- merge(result, machineDaySensorCounts, all.x = TRUE, all.y = TRUE)
    }
  }
  result[, constant := apply(.SD, 1, all),
         .SDcols = grep("constant", colnames(result), value = TRUE)]
  return(result)
}

# Finds values if a numeric vector which occur with at least a certain frequency
getBaseValues <- function(x, minFrequency = .05) {
  freqTable <- table(x) / length(x)
  return(as.numeric(names(freqTable[freqTable >= minFrequency])))
}

# Find the file in the "inputDataDir" which should contain the raw data (simply
# the biggest file) and extracts all machine-sensor combinations
summarizeMachinesAndSensorsForAggDataDir <- function(inputDataDir) {
  aggDataFiles <- list.files(path = inputDataDir, pattern = "\\.rds$", full.names = TRUE)
  filesSizes <- sapply(aggDataFiles, file.size)
  rawDataFile <- aggDataFiles[which.max(filesSizes)]
  aggDataTable <- readRDS(rawDataFile)
  return(aggDataTable[, 1, by = .(machine_name, sensor_name)][, -"V1"])
}

# Checks for intervals of values for the "attribute" which are unlike the base level
# and returns a data.table with start and end dates of such potential time series,
# requiring at least "splitThreshold" seconds between two time series candidates (if
# there is no measurement for this timespan at all, we split, too) and each time series
# having a length of "minLength" in seconds.
# If there is no base level, we simply use full days.
# "epsilon" represents the precision when checking for identity.
# The base values can be provided as parameter (if multiple ones, then base
# interval from lowest to highest creates) or the function searches for frequent
# values in the (hopefully numeric, not discrete!) data on its own.
summarizeNonBaseTimeSeriesForSingleSensor <- function(oneSensorTable, attribute = "mean",
    splitThreshold = 300, minLength = 300, epsilon = 1e-8, sensorBaseValues = NA) {
  oneSensorTable <- oneSensorTable[, c("ts", "ts_numeric", "day", attribute), with = FALSE]
  setorder(oneSensorTable, ts_numeric)
  if (is.na(sensorBaseValues)) {
    sensorBaseValues <- getBaseValues(oneSensorTable[, get(attribute)], minFrequency = .05)
  }
  # Points from the aggregation table actually represent intervals, so to get the
  # true end points (not only beginning of end interval), we consider sampling interval
  samplingInterval <- modeOfData(oneSensorTable[, diff(ts_numeric)])
  if (length(sensorBaseValues) == 0) { # default to whole days
    result <- oneSensorTable[, .(startTime = ts[1], startTime_numeric = ts_numeric[1],
        endTime = ts[.N] + samplingInterval, endTime_numeric = ts_numeric[.N] + samplingInterval),
        by = .(startDay = day)]
    setcolorder(result, c(2:5, 1)) # move startDay to position 5 (as in "ordinary"
    # result table/else-branch)
    cat("[INFO] No base value found in time series, defaulting to one-day time series.\n")
  } else {
    # Remove base value points
    if (length(sensorBaseValues) == 1) {
      oneSensorTable <- oneSensorTable[abs(get(attribute) - sensorBaseValues) > epsilon]
    } else {# build kind of base interval
      minBaseValue <- min(sensorBaseValues)
      maxBaseValue <- max(sensorBaseValues)
      oneSensorTable <- oneSensorTable[get(attribute) < minBaseValue - epsilon |
                                         get(attribute) > maxBaseValue + epsilon]
      cat(paste0("[WARN] Found more than one potential base value (",  paste0(sensorBaseValues,
          collapse = ","), "); you might want to check if the default segmentation ",
          "works for your dataset.\n"))
    }
    if (nrow(oneSensorTable) == 0) {
      cat("[INFO] One sensor does not contain any time series at all (is not ",
          "over the base value long enough).\n")
      return(data.table(startTime = as.POSIXct(character()), startTime_numeric = numeric(),
                        endTime = as.POSIXct(character()), endTime_numeric = numeric(),
                        startDay = character(), startTimeOfDay = integer()))
    }
    tsStartIdx <- c(1, oneSensorTable[, which(diff(ts_numeric) + samplingInterval >=
                                                splitThreshold) + 1])
    tsEndIdx <- c(shift(tsStartIdx, type = "lead")[-length(tsStartIdx)] - 1,
                  oneSensorTable[, .N])
    # Build result table, remove one sampling interval from start as well as end
    # (endTime = ts + sInterval - sInterval) because these are prone to contain
    # values between base and "working" level (capture transition) and therefore
    # are often below all other values in the resultinng time series
    result <- cbind(
      oneSensorTable[tsStartIdx, .(startTime = ts + samplingInterval,
                                   startTime_numeric = ts_numeric + samplingInterval)],
      oneSensorTable[tsEndIdx, .(endTime = ts, endTime_numeric = ts_numeric)]
    )
    result[, startDay := format(startTime, "%Y-%m-%d")]
  }
  result <- result[endTime_numeric - startTime_numeric >= minLength]
  result[, startTimeOfDay := 1:.N, by = startDay]
  return(result)
}

# Finds the biggest *.rds file in the "inputDataDir" and creates a summary table
# describing time series starts and ends, removing periods which are constant
# for the "attribute" and at least "splitThreshold" seconds; does only consider
# time series of at lest "minLength" seconds.
summarizeNonBaseTimeSeriesForAggDataDir <- function(inputDataDir, attribute = "mean",
    splitThreshold = 300, minLength = 300) {
  aggDataFiles <- list.files(path = inputDataDir, pattern = "\\.rds$", full.names = TRUE)
  filesSizes <- sapply(aggDataFiles, file.size)
  rawDataFile <- aggDataFiles[which.max(filesSizes)]
  aggDataTable <- readRDS(rawDataFile)
  return(aggDataTable[, summarizeNonBaseTimeSeriesForSingleSensor(oneSensorTable = .SD,
                          attribute = attribute, splitThreshold = splitThreshold,
                          minLength = minLength),
                      by = .(machine_name, sensor_name)])
}

# Returns the number of rows whose preceding row belongs to the same machine,
# same sensor and is exactly one time interval earlier (useful for forecasting)
getActualIndicesForAggDataSummary <- function(aggDataSummaryTable, intervalLength = 86400) {
  intervalsPerDay <- 86400 / intervalLength
  # Same machine, same sensor, previous interval is either one day ago (if current
  # interval is first of day) or on the same day
  return(aggDataSummaryTable[, which(
    machine_name == shift(machine_name, type = "lag") &
    sensor_name == shift(sensor_name, type = "lag") &
    ((as.numeric(difftime(day, shift(day, type = "lag"), units = "days")) == 1 &
        timeOfDay == 1 & shift(timeOfDay, type = "lag") == intervalsPerDay) |
      (day == shift(day, type = "lag") & timeOfDay == shift(timeOfDay, type = "lag") + 1))
    )])
}

# Returns the number of rows whose subsequent row belongs to the same machine,
# same sensor and is exctly one time interval later (useful for forecasting)
getPreviousIndicesForAggDataSummary <- function(aggDataSummaryTable, intervalLength = 86400) {
  intervalsPerDay <- 86400 / intervalLength
  return(aggDataSummaryTable[, which(
    machine_name == shift(machine_name, type = "lead") &
    sensor_name == shift(sensor_name, type = "lead") &
    ((as.numeric(difftime(shift(day, type = "lead"), day, units = "days")) == 1 &
        timeOfDay == intervalsPerDay & shift(timeOfDay, type = "lead") == 1) |
      (day == shift(day, type = "lead") & timeOfDay == shift(timeOfDay, type = "lead") - 1))
  )])
}

# Only retains intervals in the table with the previous or next time interval
# of the same machine and sensor also existing
removeIntervalsWithoutPredAndSucc <- function(aggDataSummaryTable,
                                              intervalLength = 86400) {
  result <- copy(aggDataSummaryTable)
  setorder(result, machine_name, sensor_name, day, timeOfDay)
  return(result[sort(union(getActualIndicesForAggDataSummary(result, intervalLength),
                      getPreviousIndicesForAggDataSummary(result, intervalLength)))])
}

# Analyzes a table containing machine_name, sensor_name, day, time of day for
# consecutive time intervals of the same machine and sensor, returning the
# corresponding indices with their predecessors as list (which is required for
# forecasting -> see ClusteringEvaluationUtility)
createPrevActListForAggDataSummary <- function(aggDataSummaryTable,
                                               intervalLength = 86400) {
  result <- list(prev = getPreviousIndicesForAggDataSummary(aggDataSummaryTable,
                    intervalLength = intervalLength),
                 actual = getActualIndicesForAggDataSummary(aggDataSummaryTable,
                    intervalLength = intervalLength))
  if (any(result$act != result$prev + 1)) {
    stop("Error in computation of prevActList")
  }
  return(result)
}

# Reads all data.tables representing Spark aggregations as RDS file from the
# "inputDataDir" and converts the "aggColumnList" (each element being a character
# vector, resulting in univariate time series if one element and multivariate if
# several elements; default: all columns converted to univariate time series) to
# a list of time series which is saved in the "outputDataDir";
# - "machineDaySensorCombinations" can be used to limit data, e.g. certain machines,
# sensors days etc.
# - "timeSeriesCombinations" can be used to create variable length time series based
# on a table with machine/sensor names, start and end dates
# - "from" and "to" specify start and end time
# - if "rawAggColumn" is not NA, then only one attribute will be used for the
# biggest dataset and saved as "Raw"
# - "intervalLength" determines if days should be split into multiple time series,
# each having this length (in seconds; default: one day is one time series)
saveAllAggregateTimeSeries <- function(datasetName, inputDataDir = "./",
    outputDataDir = "./aggregations/", machineDaySensorCombinations = NULL,
    timeSeriesCombinations = NULL, intervalLength = 86400, aggColumnList = NULL,
    rawAggColumn = NA, naReplace = NA, normalize = "") {
  if (!is.null(machineDaySensorCombinations) && nrow(machineDaySensorCombinations) == 0) {
    stop(paste("Cannot create time series.",
               "The passed table of machine-sensor-day combinations is empty."))
  }
  if (!is.null(timeSeriesCombinations) && nrow(timeSeriesCombinations) == 0) {
    stop(paste("Cannot create time series.",
               "The passed table of time series combinations is empty."))
  }
  dir.create(outputDataDir, showWarnings = FALSE, recursive = TRUE)
  tsListNames <- NULL
  if (!is.na(rawAggColumn)) { # determine file with most data
    aggFileSizes <- sapply(list.files(path = inputDataDir, pattern = "\\.rds$",
                                      full.names = TRUE), file.size)
    rawDataFile <- basename(names(aggFileSizes)[which.max(aggFileSizes)])
  }
  for (fileName in list.files(path = inputDataDir, pattern = "rds$")) {
    cat("##### File: ", fileName, "######\n")
    aggDataTable <- readRDS(file = paste0(inputDataDir, fileName))
    if (intervalLength > 0) {
      # addTimeOfDay() only addresses time series of fixed length
      aggDataTable <- addTimeOfDay(aggDataTable, intervalLength = intervalLength)
    }
    if (!is.na(rawAggColumn) && fileName == rawDataFile) {
      setnames(aggDataTable, rawAggColumn, "Raw")
      currAggColumnList <- list("Raw")
    } else if (is.null(aggColumnList)) {
      currAggColumnList <- as.list(setdiff(colnames(aggDataTable),
          c("machine_name", "sensor_name", "ts", "ts_numeric", "day", "secOfDay", "timeOfDay")))
    } else {
      currAggColumnList <- aggColumnList
    }
    if (!is.null(machineDaySensorCombinations)) {
      aggDataTable <- merge(machineDaySensorCombinations, aggDataTable, all.x = TRUE)
    }
    if (!is.na(naReplace)) {
      aggDataTable[is.na(aggDataTable)] <- naReplace
    }
    aggDataTable <- interpolateAllAggData(aggDataTable, onlyExistingCombinations = TRUE,
                                          intervalLength = intervalLength)
    samplingInterval <- modeOfData(aggDataTable[, diff(sort(unique(ts_numeric)))])
    for (aggColumns in currAggColumnList) {
      # Aggregation name should not contain underscore (has special meaning in other scripts)
      aggName <- paste(aggColumns, collapse = "-")
      cat(aggName, "|", sep = "")
      if (intervalLength > 0) {
        tsList <- convertToFixedIntervalTimeSeriesList(aggDataTable,
            attributeNames = aggColumns, normalize = normalize)
        tsLength <- unique(sapply(tsList, function(ts) {
          return(ifelse(is.matrix(ts), nrow(ts), length(ts)))
        }))
        if (length(tsLength) > 1) {
          stop("Error - time series of different length were created.")
        }
      } else {
        tsList <- convertToVariableLengthTimeSeriesList(aggDataTable,
            timeSeriesCombinations = timeSeriesCombinations,
            attributeNames = aggColumns, normalize = normalize)
        # Time series have variable length, but for evaluation purpose we need
        # a level, so we use the number of points that would be in a series
        # representing one day
        tsLength <- 86400 / samplingInterval
      }
      if (is.null(tsListNames)) {
        tsListNames <- names(tsList)
      } else if ((length(tsList) != length(tsListNames)) ||
                 !all(names(tsList) == tsListNames)) {
        stop("Number and/or names of time series have changed compared to first
             aggregation level/type.")
      }
      saveRDS(tsList, file = paste0(outputDataDir, datasetName, "_", aggName,
                                    "_", tsLength, ".rds"))
    }
    cat("\n")
  }
  return(TRUE)
}

### Index vector methods for aggregated data time series ###

# Dispatching routine to routine a time series list which can also read files
# or search directories (just for convenience, error handling could be improved)
retrieveTSList <- function(tsListOrPath) {
  if (is.list(tsListOrPath)) {
    return(tsListOrPath)
  } else if (file.exists(tsListOrPath)) {
    tsList <- readRDS(tsListOrPath)
    if (is.list(tsList)) {
      return(tsList)
    } else {
      stop("Found file, but is not a list.")
    }
  } else if (dir.exists(tsListOrPath)) {
    fileName <- list.files(tsListOrPath, pattern = "\\.rds", full.names = TRUE)[1]
    tsList <- readRDS(fileName)
    if (is.list(tsList)) {
      return(tsList)
    } else {
      stop("Used first *.rds file in directory, but is not a list.")
    }
  } else {
    stop("Parameter is neither a list nor a file nor a directory.")
  }
}

# Retrieves a certain part from names of a list (as created by the
# convertToTimeSeriesList.aggData method) and converts it to an integer vector
getIndexVectorFromTSList <- function(tsListOrPath, namePart = 1, sep = "-") {
  tsList <- retrieveTSList(tsListOrPath)
  relevantNamePart <- sapply(strsplit(names(tsList), split = sep),
                             function(x) x[[namePart]])
  return(as.integer(as.factor(relevantNamePart)))
}

# Creates an integer vector marking different machines in a list of time series
# based on its names (convertToTimeSeriesList.aggData method).
getMachineIndexVector <- function(tsListOrPath) {
  getIndexVectorFromTSList(tsListOrPath, namePart = 1)
}

# Creates an integer vector marking different sensors in a list of time series
# based on its names (convertToTimeSeriesList.aggData method).
getSensorIndexVector <- function(tsListOrPath) {
  getIndexVectorFromTSList(tsListOrPath, namePart = 2)
}

# Creates an integer vector marking different days in a list of time series
# based on its names (convertToTimeSeriesList.aggData method).
getDayIndexVector <- function(tsListOrPath) {
  getIndexVectorFromTSList(tsListOrPath, namePart = 3)
}

# Creates an integer vector marking different intervals of days in a list of
# time series based on its names (convertToTimeSeriesList.aggData method).
getTimeOfDayIndexVector <- function(tsListOrPath) {
  getIndexVectorFromTSList(tsListOrPath, namePart = 4)
}

# Creates an integer vector marking different days of the week in a list of
# time series based on its names (convertToTimeSeriesList.aggData method).
getDayOfWeekIndexVector <- function(tsListOrPath, sep = "-") {
  tsList <- retrieveTSList(tsListOrPath)
  dateString <- sapply(strsplit(names(tsList), split = sep),
                       function(x) x[[3]])
  dayOfWeek <- format(as.POSIXct(dateString, format = "%Y_%m_%d", tz = "GMT"), "%u")
  return(as.integer(dayOfWeek))
}

# Creates an integer vector marking weekday (1) vs weekend (2) in a list of
# time series based on its names (convertToTimeSeriesList.aggData method).
getWeekDayVsEndIndexVector <- function(tsListOrPath, sep = "-") {
  dayOfWeek <- getDayOfWeekIndexVector(tsListOrPath, sep = sep)
  result <- dayOfWeek %in% 6:7 + 1
  return(as.integer(result))
}
