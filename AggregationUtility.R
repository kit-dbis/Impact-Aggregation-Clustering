library(batchtools)
library(data.table)
library(FastTSDistances)
library(snow) # necessary for batchtools, but not imported there

### Functions ###

# Calculates aggregation levels for piecewise aggregation using powers of 2 from
# minRepresentationLength (be aware that some distance measures might require a
# minimum time series length greater than 1) to a maximum length (which depends
# on the length of the shortest time series and a minimum length of the aggregation
# window)
calcAggLevels.PAA <- function(tsList, minRepresentationLength = 4,
    minWindowLength = 1) {
  shortestTSLength <- min(sapply(tsList, function(x) {
    return(ifelse(is.matrix(x), nrow(x), length(x)))
  }))
  maxRepresentationLength <- floor(shortestTSLength / minWindowLength)
  aggLevels <- 2^(floor(log(maxRepresentationLength, 2)):ceiling(log(minRepresentationLength, 2)))
  return(aggLevels)
}

# Calculates piecewise aggregation by dividing each of the original series
# into "aggLevel" equal-length segments which are represented by one aggregate
# each; if "baseLength" is provided, this indicates variable length time series
# and the concrete aggregation level will be computed as the "aggLevel" value
# multiplied with the ratio of actual length and base length
calcAggregation.PAA <- function(tsList, aggLevel, baseLength = NA,
                                FUN = PAA_fast, ...) {
  additionalParams <- list(...)
  return(lapply(tsList, function(ts) {
    if (is.matrix(ts)) {
      if (!is.na(baseLength)) {
        windowCount <- ceiling(nrow(ts) / baseLength * aggLevel)
      } else {
        windowCount <- aggLevel
      }
      result <- cbind(sapply(1:ncol(ts), function(j) {
        do.call(FUN, args = c(list(ts[, j]), windowCount = windowCount, additionalParams))
      }))
      colnames(result) <- colnames(ts)
      return(result)
    } else {
      if (!is.na(baseLength)) {
        windowCount <- ceiling(length(ts) / baseLength * aggLevel)
      } else {
        windowCount <- aggLevel
      }
      return(do.call(FUN, args = c(list(ts), windowCount = windowCount, additionalParams)))
    }
  }))
}

# Dummy for raw data (no aggregation to compute)
calcAggregation.raw <- function(tsList, aggLevel, baseLength = NA) {
  return(tsList)
}

# Dummy, there are no real levels for raw data
calcAggLevels.raw <- function(tsList) {
  tsLengths <- sapply(tsList, function(ts) {
    return(ifelse(is.matrix(ts), nrow(ts), length(ts)))
  })
  if (any(tsLengths != tsLengths[1])) {
    warning(paste0("Time series do not have the same length, we take the ",
                    "maximum as raw aggregation level."))
  }
  return(max(tsLengths))
}

### Aggregations with pre-configured parameters ###

DEFAULT_AGGREGATIONS <- list(
  Raw = list(
    name = "Raw", levelMethod = calcAggLevels.raw, levelParams = list(),
    aggMethod = calcAggregation.raw, aggParams = list()
  ),
  PAA.Mean = list(
    name = "PAA.Mean", levelMethod = calcAggLevels.PAA, levelParams = list(),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PAA_fast)
  ),
  PAA.Median = list(
    name = "PAA.Median", levelMethod = calcAggLevels.PAA, levelParams = list(),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PMedAA_fast)
  ),
  PAA.Min = list(
    name = "PAA.Min", levelMethod = calcAggLevels.PAA, levelParams = list(),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PMinAA_fast)
  ),
  PAA.Max = list(
    name = "PAA.Max", levelMethod = calcAggLevels.PAA, levelParams = list(),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PMaxAA_fast)
  ),
  PAA.SD = list(
    name = "PAA.SD", levelMethod = calcAggLevels.PAA, levelParams = list(),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PSDAA_fast, sample = TRUE)
  ),
  PAA.Skew = list(
    name = "PAA.Skew", levelMethod = calcAggLevels.PAA, levelParams = list(minWindowLength = 2),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PSkewAA_fast, nanReplace = 0)
  ),
  PAA.Kurt = list(
    name = "PAA.Kurt", levelMethod = calcAggLevels.PAA, levelParams = list(minWindowLength = 2),
    aggMethod = calcAggregation.PAA, aggParams = list(FUN = PKurtAA_fast, nanReplace = 0, excess = TRUE)
  )
)

# Takes aggregation names as used in Spark (plus the "Raw" level which we might
# add manually) and converts them to the names in our default list of aggregations
sparkAggNamesToLocalAggNames <- function(sparkAggNames) {
  return(sapply(sparkAggNames, switch, Raw = "Raw", mean = "PAA.Mean",
                median = "PAA.Median", min = "PAA.Min", max = "PAA.Max",
                stddev = "PAA.SD", skewness = "PAA.Skew", kurtosis = "PAA.Kurt", ""))
}

# Calculates multiple aggregations from "aggNames" for a list of time series,
# returning a matrix for each object from "tsList"; the input matrices can be
# multivariate, so each aggregation will be computed for each attribute
calcMultiAggregation.PAA <- function(tsList, aggLevel, baseLength = NA, aggNames) {
  # Perform each aggregation independently (time series can be uni- or multivariate,
  # does not matter at this point [aggMethod should be able to handle both])
  univariateAggregations <- lapply(aggNames, function(aggName) {
    do.call(DEFAULT_AGGREGATIONS[[aggName]]$aggMethod, c(tsList = list(tsList),
        aggLevel = aggLevel, baseLength = baseLength, DEFAULT_AGGREGATIONS[[aggName]]$aggParams))
  })
  names(univariateAggregations) <- sapply(DEFAULT_AGGREGATIONS[aggNames],
                                          function(x) x$name)
  if (is.matrix(tsList[[1]])) { #multivariate time series
    aggMatrixColumnNames <- paste0(rep(colnames(tsList[[1]]), times = length(aggNames)),
        ".", rep(names(univariateAggregations), each = ncol(tsList[[1]])))
    # Combine the different aggregation types for each time series to a matrix
    multivariateAggregations <- lapply(1:length(tsList), function(i) {
      aggMatrix <- do.call(cbind, lapply(names(univariateAggregations), function(aggName)
        univariateAggregations[[aggName]][[i]]))
      colnames(aggMatrix) <- aggMatrixColumnNames
      return(aggMatrix)
    })
  } else {#univariate time series
    # Combine the different aggregation types for each time series to a matrix
    multivariateAggregations <- lapply(1:length(tsList), function(i)
      cbind(sapply(names(univariateAggregations), function(aggName)
        univariateAggregations[[aggName]][[i]])))
  }
  names(multivariateAggregations) <- names(tsList)
  return(multivariateAggregations)
}

# Takes a list of character vectors representing Spark aggregation types and a
# vector of aggregation levels to return functions for local aggregations. All
# aggregations (except "Raw") get the same aggregation levels. Using multivariate
# aggregations (like mean+min+max) is possible
# "variableLength" indicates that time series of different lengths should be
# aggregated which requires a more flexible handling of aggregation levels.
sparkAggNameListToAggFunctionList <- function(sparkAggNameList, aggregationLevels,
                                              variableLength = TRUE) {
  # We need to force evaluation as the aggregation levels might be created by a
  # function which might not be available later (e.g. in a batchtools process)
  aggregationLevels <- force(aggregationLevels)
  if (variableLength) {
    baseLength <- 2880 # 30 secs
  } else {
    baseLength <- NA
  }
  aggFunctionList <- lapply(sparkAggNameList, function(sparkAggNames) {
    localAggNames <- sparkAggNamesToLocalAggNames(sparkAggNames)
    if (any(localAggNames == "")) {
      cat("[WARN] The aggregate(s) ", paste0("\"", sparkAggNames[localAggNames == ""], "\"", collapse = ", "),
          " cannot be computed locally. If they are contained in your dataset, you can ignore this warning.\n", sep = "")
    }
    localAggNames <- localAggNames[localAggNames != ""]
    if (length(localAggNames) == 0) { # no local aggregate computation possible
      return(NULL)
    } else if (length(localAggNames) == 1) { # univariate aggregation (simple case)
      result <- DEFAULT_AGGREGATIONS[[localAggNames]]
      if (localAggNames != "Raw") {
        result$levelMethod = function(tsList, ...) return(aggregationLevels)
        result$aggParams[["baseLength"]] <- baseLength
      }
    } else {# multivariate aggregation (we need to combine function calls)
      result <- list(
        name = paste0("PAA.", paste0(gsub("PAA.", "", localAggNames, fixed = TRUE), collapse = "-")),
        levelMethod = function(tsList, ...) return(aggregationLevels),
        levelParams = list(),
        aggMethod = calcMultiAggregation.PAA,
        aggParams = list(baseLength = baseLength, aggNames = localAggNames)
      )
    }
    return(result)
  })
  names(aggFunctionList) <- sapply(sparkAggNameList, paste0, collapse = "-")
  aggFunctionList <- aggFunctionList[!sapply(aggFunctionList, is.null)]
  return(aggFunctionList)
}

### Whole dataset aggregation ###

# Stops if assumptions of the aggregation routine are violated. It is important
# that list elements are named such that datasets and aggregations can be identified
# correctly.
checkAggregationAssumptions <- function(tsListList, aggList, outputBasePathList) {
  if (is.null(names(tsListList)) ||
      any(sapply(names(tsListList), is.null)) ||
      any(sapply(names(tsListList), is.na)) ||
      any(names(tsListList) == "")) {
    stop("All your datasets should be named.")
  }
  if (is.null(names(aggList)) ||
      any(sapply(names(aggList), is.null)) ||
      any(sapply(names(aggList), is.na)) ||
      any(names(aggList)  == "")) {
    stop("All your aggregations should be named.")
  }
  if (length(outputBasePathList) > 0 &&
      (is.null(names(outputBasePathList)) ||
       any(sapply(names(outputBasePathList), is.null)) ||
       any(sapply(names(outputBasePathList), is.na)) ||
       any(names(outputBasePathList) == "") ||
       length(setdiff(names(outputBasePathList), names(tsListList))) > 0)) {
    stop("All your output directory paths should be named like the
         corresponding dataset in tsListList.")
  }
}

# Computes all aggregations from the "aggList" for all lists of time series in the
# "tsListList", choosing appropriate levels (depending on the dataset and aggregation
# method). You can provide a named list containing an output directory for each
# dataset (element of "tsListList"), otherwise the outputs are saved in the working
# directory. All aggregations are saved as .rds file, the function always returns
# TRUE.
aggregateAndSaveData <- function(tsListList, aggList, outputBasePathList = list(),
    expRegistry = NULL, cpuCores = 8) {
  checkAggregationAssumptions(tsListList, aggList, outputBasePathList)
  if (is.null(expRegistry)) {
    expRegistry <- makeExperimentRegistry(file.dir = NA, seed = 1)
  }
  expRegistry$cluster.functions <- makeClusterFunctionsSocket(ncpus = cpuCores)
  expRegistry$source <- "AggregationUtility.R" # file will be sourced on start
  # Problems: different datasets (which might have individual output paths)
  for (i in 1:length(tsListList)) {
    datasetName <- names(tsListList)[i]
    outputBasePath <- outputBasePathList[[datasetName]]
    if (is.null(outputBasePath)) {
      outputBasePath <- "./"
    }
    if (!dir.exists(outputBasePath)) {
      dir.create(outputBasePath, recursive = TRUE)
    }
    addProblem(
      name = datasetName,
      data = list(tsList = tsListList[[i]], aggList = aggList, outputBasePath = outputBasePath),
      fun = NULL, # no pre-processing needed, aggregation is part of experiment
      reg = expRegistry
    )
  }
  # Algorithms: different aggregations and levels
  addAlgorithm(name = "Aggregation", fun = function(data, job, instance, ...) {
    cat("Aggregating with ", job$pars$algo.pars$aggName, "...\n", sep = "")
    aggName <- job$pars$algo.pars$aggName
    aggLevels <- do.call(
      what = data$aggList[[aggName]]$levelMethod,
      args = c(tsList = list(data$tsList), aggList[[aggName]]$levelParams)
    )
    for (aggLevel in aggLevels) {
      cat("Working on level ", aggLevel, "...\n", sep = "")
      aggTSList <- do.call(
        what = data$aggList[[aggName]]$aggMethod,
        args = c(tsList = list(data$tsList), aggLevel = aggLevel, data$aggList[[aggName]]$aggParams)
      )
      saveRDS(aggTSList, file = paste0(data$outputBasePath, job$prob.name, "_",
                                       aggName, "_", aggLevel, ".rds"))
    }
    return(TRUE)
  }, reg = expRegistry)
  # Experiment: Combine problems + algorithms
  addExperiments(
    prob.designs = NULL,
    algo.designs = list(Aggregation = data.table(aggName = names(aggList))),
    reg = expRegistry
  )
  submitJobs(reg = expRegistry)
  waitForJobs(reg = expRegistry)
  return(TRUE)
}
