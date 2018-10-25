library(corrplot)

source("BatchUtility.R")
source("GeneralPlotRoutinesUtility.R")

# Retrieves the biggest dataset from a batch experiment
readBaseDataset <- function(dataPath) {
  biggestTSFile <- names(which.max(sapply(
    list.files(paste0(dataPath, "aggData/"), full.names = T), file.size)))
  return(readRDS(biggestTSFile))
}

# Plots the time series at the lowest aggregation (aka base level)
saveExperimentBaseTimeSeriesPlots <- function(dataPath) {
  biggestTSFile <- names(which.max(sapply(
    list.files(paste0(dataPath, "aggData/"), full.names = T), file.size)))
  tsList <- readRDS(biggestTSFile)
  # Retain only "Aggregation_Level" (not dataset, not file ending)
  outputDir <- paste0(gsub("\\.rds", "", gsub("^[^_]+_", "", basename(biggestTSFile))), "/")
  plotTimeSeriesList(tsList, paste0(gsub("data", "plots", dataPath), outputDir))
  return(TRUE)
}

# Plots the original time series, copies them according to the clustering and creates
# cluster overview plots for a certain aggregation-level-dissimilarity-clustering
# combination (you can also pass character vectors longer than one which results
# in the cross-product to be plotted). Necessary paths are built dynamically from
# the "dataPath" and old plots are deleted.
saveExperimentClusterPlots <- function(dataPath, clusterAssignmentsTable = NULL, aggregation,
    level, dissimilarity, clustering) {
  if (is.null(clusterAssignmentsTable)) {
    clusterAssignmentsTable <- readRDS(paste0(dataPath, "results/clusteringResult.rds"))
  }
  # Cross-product
  settingsTable <- data.table(expand.grid(aggregation, level, dissimilarity, clustering,
                                          stringsAsFactors = FALSE))
  if (nrow(settingsTable) > 1) {
    for (i in 1:nrow(settingsTable)) {
      cat("--Working on ", paste(settingsTable[i,], collapse = "-"), " (combination ",
          i, "/", nrow(settingsTable), ")--\n", sep = "")
      saveExperimentClusterPlots(dataPath = dataPath, clusterAssignmentsTable = clusterAssignmentsTable,
        aggregation = settingsTable[i, Var1], level = settingsTable[i, Var2],
        dissimilarity = settingsTable[i, Var3], clustering = settingsTable[i, Var4])
    }
    return(TRUE)
  }
  plotPath <- gsub("data", "plots", dataPath)
  relevantAssignmentsTable <- clusterAssignmentsTable[Aggregation == aggregation &
      Level == level & Dissimilarity == dissimilarity & Clustering == clustering]
  if (nrow(relevantAssignmentsTable) != 1) {
    stop(paste0("Your combination of aggregation-level-dissimilarity-clustering ",
                "does not exist or is not unique."))
  }
  dataset <- relevantAssignmentsTable[, Dataset]
  assignments <- relevantAssignmentsTable[,
      as.integer(mget(grep("A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE)))]
  tsList <- readRDS(list.files(path = paste0(dataPath, "aggData/"),
      pattern = paste0(aggregation, "_", level, "\\.rds"), full.names = TRUE))
  originalTSPlotPath <- paste0(plotPath, aggregation, "_", level, "/")
  existingTsPlotFiles <- list.files(originalTSPlotPath, pattern = "\\.png", full.names = TRUE)
  if (length(existingTsPlotFiles) != length(tsList) ||
      !all(basename(existingTsPlotFiles) == paste0(names(tsList), ".png"))) {
    if (length(existingTsPlotFiles) != 0) {
      cat("Deleting previous time series plots as their names did not match the",
          "time series list.\n")
      unlink(existingTsPlotFiles)
    }
    cat("Plotting unclustered time series.\n")
    plotTimeSeriesList(tsList, directory = originalTSPlotPath)
  } else {
    cat("Using existing plots of original time series.\n")
  }
  clusterPlotPath <- paste0(originalTSPlotPath, dissimilarity, "_", clustering)
  if (length(list.files(clusterPlotPath, include.dirs = TRUE)) > 0) {
    cat("Deleting previous cluster plots.\n")
    unlink(clusterPlotPath, recursive = TRUE) # delete old
  }
  clusterPlotPath <- paste0(clusterPlotPath, "/") # unlink() does not work with slash at end
  cat("Creating cluster plots and copying time series plots according to clustering.\n")
  if (clustering == "AffProp") {
    newAssignments <- clusterTSWithAffinityPropagation(tsList = tsList,
        dissMatrix = readRDS(paste0(dataPath, "dissMatrices/DissMatrix_", dataset, "_",
                                    dissimilarity, "_", aggregation, "_", level, ".rds")),
        tsPlotPath = originalTSPlotPath, clusPlotPath = clusterPlotPath)
    if (any(newAssignments != assignments)) {
      cat("Tried to re-run clustering algorithm, but got different results.\n")
    }
  } else if (clustering == "PAM") { # re-run to get centroids
    newAssignments <- clusterTSWithPAM(tsList = tsList,
        dissMatrix = readRDS(paste0(dataPath, "dissMatrices/DissMatrix_", dataset, "_",
                                    dissimilarity, "_", aggregation, "_", level, ".rds")),
        tsPlotPath = originalTSPlotPath, clusPlotPath = clusterPlotPath,
        k = DEFAULT_CLUSTERINGS[["PAM"]]$params$k)
    if (any(newAssignments != assignments)) { # re-run to get centroids
      cat("Tried to re-run clustering algorithm, but got different results.\n")
    }
  } else {# routines without centroids
    saveClusterPlots(tsList = tsList, assignments = assignments, tsPlotPath = originalTSPlotPath,
                     clusPlotPath = clusterPlotPath)
  }
  return(TRUE)
}

# Creates forecast plots for a certain aggregation-level-dissimilarity-clustering
# combination (you can also pass a vector for each of the parameters and the cross-
# product will be plotted) and saves them on the hard drive. Base (with the raw
# time series) or relative (with the aggregated time series) forecasting possible.
saveExperimentForecastPlots <- function(dataPath, aggregation, level, dissimilarity,
    clustering, relative = TRUE, withBackground = TRUE) {
  # Cross-product
  settingsTable <- data.table(expand.grid(aggregation, level, dissimilarity, clustering,
      relative, stringsAsFactors = FALSE))
  if (nrow(settingsTable) > 1) {
    return(rbindlist(lapply(seq_len(nrow(settingsTable)), function(i) {
      cat("--Working on ", paste(settingsTable[i,], collapse = "-"), " (combination ",
          i, "/", nrow(settingsTable), ")--\n", sep = "")
      saveExperimentForecastPlots(dataPath = dataPath, aggregation = settingsTable[i, Var1],
          level = settingsTable[i, Var2], dissimilarity = settingsTable[i, Var3],
          clustering = settingsTable[i, Var4], relative = settingsTable[i, Var5],
          withBackground = withBackground)
    })))
  }
  cat("Loading data.\n")
  trainAssignmentsTable <- readRDS(paste0(dataPath, "results/clusteringResult.rds"))
  relevantTrainAssignmentsTable <- trainAssignmentsTable[Aggregation == aggregation &
      Level == level & Dissimilarity == dissimilarity & Clustering == clustering]
  if (nrow(relevantTrainAssignmentsTable) != 1) {
    stop(paste0("Your combination of aggregation-level-dissimilarity-clustering ",
                "does not exist or is not unique for the training data."))
  }
  dataset <- relevantTrainAssignmentsTable[, Dataset]
  testAssignmentsTable <- readRDS(paste0(dataPath, "results/testDataAssignmentResult.rds"))
  relevantTestAssignmentsTable <- testAssignmentsTable[Aggregation == aggregation &
      Level == level & Dissimilarity == dissimilarity & Clustering == clustering]
  if (nrow(relevantTestAssignmentsTable) != 1) {
    stop(paste0("Your combination of aggregation-level-dissimilarity-clustering ",
                "does not exist or is not unique for the test data."))
  }
  trainDataPrevActList <- readRDS(paste0(dataPath, "results/trainDataPrevActList.rds"))
  testDataPrevActList <- readRDS(paste0(dataPath, "results/testDataPrevActList.rds"))
  if (relative) { # we make forecasts always on "mean" aggregation (see BatchUtility.R)
    trainDataset <- readRDS(paste0(dataPath, "fullTrainAggData/", dataset, "_mean_",
                                   level, ".rds"))
    testDataset <- readRDS(paste0(dataPath, "fullTestAggData/", dataset, "_mean_",
                                  level, ".rds"))
  } else {
    rawTrainDatasetFiles <- list.files(path = paste0(dataPath, "fullTrainAggData/"),
        pattern = paste0(dataset, "_Raw_[0-9]+\\.rds"), full.names = TRUE)
    if (length(rawTrainDatasetFiles) == 1) {
      trainDataset <- readRDS(rawTrainDatasetFiles)
    } else {
      stop("Found multiple base (raw) train datasets.")
    }
    rawTestDatasetFiles <- list.files(path = paste0(dataPath, "fullTestAggData/"),
        pattern = paste0(dataset, "_Raw_[0-9]+\\.rds"), full.names = TRUE)
    if (length(rawTestDatasetFiles) == 1) {
      testDataset <- readRDS(rawTestDatasetFiles)
    } else {
      stop("Found multiple base (raw)test datasets.")
    }
  }
  forecastPlotDir <- paste0(gsub("data", "plots", dataPath), "Forecasts/", aggregation, "_",
      level, "_", dissimilarity, "_", clustering, "_", ifelse(relative, "relative", "base"))
  if (length(list.files(forecastPlotDir, include.dirs = TRUE)) > 0) {
    cat("Deleting previous forecast plots.\n")
    unlink(forecastPlotDir, recursive = TRUE)
  }
  forecastPlotDir <- paste0(forecastPlotDir, "/")
  cat("Creating forecast plots.\n")
  plotForecasts(
    newDatasetPrev = testDataset[testDataPrevActList[[1]]$prev],
    newDatasetActual = testDataset[testDataPrevActList[[1]]$actual],
    oldDatasetPrev = trainDataset[trainDataPrevActList[[1]]$prev],
    oldDatasetActual = trainDataset[trainDataPrevActList[[1]]$actual],
    newDataAssignments = relevantTestAssignmentsTable[,
        as.integer(mget(grep("^A[0-9]+", colnames(relevantTestAssignmentsTable), value = TRUE)))],
    oldDataAssignments = relevantTrainAssignmentsTable[,
        as.integer(mget(grep("^A[0-9]+", colnames(relevantTrainAssignmentsTable), value = TRUE)))],
    dissMatrix = readRDS(paste0(dataPath, "testDissMatrices/DissMatrix_", dataset, "_",
                                dissimilarity, "_", aggregation, "_", level, ".rds")),
    plotDir = forecastPlotDir, dissWeightConversion = "gaussian",
    withBackground = withBackground, progressBar = FALSE
  )
  cat("Calculating forecast errors.\n")
  result <- calcForecastError(
    newDatasetPrev = testDataset[testDataPrevActList[[1]]$prev],
    newDatasetActual = testDataset[testDataPrevActList[[1]]$actual],
    oldDatasetPrev = trainDataset[trainDataPrevActList[[1]]$prev],
    oldDatasetActual = trainDataset[trainDataPrevActList[[1]]$actual],
    newDataAssignments = relevantTestAssignmentsTable[,
        as.integer(mget(grep("^A[0-9]+", colnames(relevantTestAssignmentsTable), value = TRUE)))],
    oldDataAssignments = relevantTrainAssignmentsTable[,
        as.integer(mget(grep("^A[0-9]+", colnames(relevantTrainAssignmentsTable), value = TRUE)))],
    dissMatrix = readRDS(paste0(dataPath, "testDissMatrices/DissMatrix_", dataset, "_",
        dissimilarity, "_", aggregation, "_", level, ".rds")),
    dissWeightConversion = "gaussian"
  )
  result[, Dataset := dataset]
  result[, Dissimilarity := dissimilarity]
  result[, Clustering := clustering]
  result[, Aggregation := aggregation]
  result[, Level := level]
  result[, Relative := relative]
  setcolorder(result, c(4:8, 1:3, 9))
  return(result)
}

# Prints tables of the main experiment categories for a summary table (which
# should be filtered for a criterion or the categories are simply balanced)
tableExperimentCategories <- function(clusterSummaryTable, hint = "") {
  cat("### Tables for categories")
  if (hint != "") cat(" (", hint, ")", sep = "")
  cat("\n\n")
  if (clusterSummaryTable[, uniqueN(Dataset)] > 1) {
    print(clusterSummaryTable[, sort(table(Dataset), decreasing = TRUE)])
    cat("\n")
  }
  print(clusterSummaryTable[, sort(table(Dissimilarity), decreasing = TRUE)])
  cat("\n")
  print(clusterSummaryTable[, sort(table(Clustering), decreasing = TRUE)])
  cat("\n")
  print(clusterSummaryTable[, sort(table(Aggregation), decreasing = TRUE)])
  cat("\n")
  print(clusterSummaryTable[, table(Level)])
  cat("\n")
  return(TRUE)
}

# Prints an aggregate like mean, median for any evaluation table and the main
# experiment categories
summarizeExperimentCategories <- function(clusterSummaryTable, funcName = "mean", hint = "") {
  cat("### \"", funcName, "()\" for categories", sep = "")
  if (hint != "") cat(" (", hint, ")", sep = "")
  cat("\n\n")
  if (clusterSummaryTable[, uniqueN(Dataset)] > 1) {
    print(clusterSummaryTable[, do.call(funcName, list(Value, na.rm = TRUE)), by = "Dataset"][order(-V1)])
    cat("\n")
  }
  print(clusterSummaryTable[, do.call(funcName, list(Value, na.rm = TRUE)), by = "Dissimilarity"][order(-V1)])
  cat("\n")
  print(clusterSummaryTable[, do.call(funcName, list(Value, na.rm = TRUE)), by = "Clustering"][order(-V1)])
  cat("\n")
  print(clusterSummaryTable[, do.call(funcName, list(Value, na.rm = TRUE)), by = "Aggregation"][order(-V1)])
  cat("\n")
  print(clusterSummaryTable[, do.call(funcName, list(Value, na.rm = TRUE)), by = "Level"][order(-Level)])
  cat("\n")
  return(TRUE)
}

# Prints the best experimental settings per category, keeping the other settings
# fixed; if "group" provided, best experimental settings per group and category,
# else a ranking of all experimental settings per category
summarizeBestExperimentSettings <- function(clusterSummaryTable, group = NULL,
      maximize = TRUE) {
  for (categoryName in setdiff(colnames(clusterSummaryTable), c("Value", group))) {
    numberOfSettings <- clusterSummaryTable[, uniqueN(get(categoryName))]
    if (numberOfSettings > 1) {
      cat("Category: ", categoryName, " (", numberOfSettings, " settings)\n", sep = "")
      if (is.null(group)) { # show ranking of all category values
        print(countBestExperimentSetting(clusterSummaryTable, categoryName = categoryName,
            maximize = maximize, relativeCount = T, na.rm = T))
      } else {# pick only best category value per group
        print(clusterSummaryTable[, countBestExperimentSetting(.SD,
            categoryName = categoryName, maximize = maximize, na.rm = T,
            relativeCount = T)[1], by = group])
      }
      cat("\n")
    }
  }
  return(TRUE)
}

# Calculates the correlation matrix between all values of a category (indices,
# dissimilarities etc.) and optionally creates plots
# - if the values strongly depend on a second category, you can use "byCategory":
# all correlations are calculated for each value of this attribute separately and
# then averaged
# - "infStrategy" determines how to handle Inf values, options are "drop", "max"
# (replace Inf with highest non-Inf value of that index) or doing nothing (which
# might result in Na correlations)
# - "correlationMethod" is passed to cor()
summarizeExperimentSettingCorrelation <- function(clusterSummaryTable, category,
    byCategory = NA, infStrategy = "", plot = TRUE, plotFilePath = NA, correlationMethod = "pearson") {
  clusterSummaryTable <- copy(clusterSummaryTable)
  descriptiveColumns <- setdiff(colnames(clusterSummaryTable), c("Value", category))
  clusterSummaryTable[, Setting := .GRP, by = descriptiveColumns]
  if (infStrategy == "drop") {
    clusterSummaryTable[is.infinite(Value), Value := NA]
  } else if (infStrategy == "max") {
    clusterSummaryTable[, Value := pmin(Value,
        max(Value[!is.infinite(Value)], na.rm = TRUE)), by = Index]
  }
  if (is.na(byCategory)) {
    categoryValueTable <- dcast(clusterSummaryTable, value.var = "Value",
                                formula = paste0("Setting ~ ", category))
    categoryValueTable[, Setting := NULL]
    # Ignore warning if standard deviation is zero
    corrMatrix <- suppressWarnings(cor(categoryValueTable, use = "pairwise.complete.obs",
                                       method = correlationMethod))
  } else {# correlate for each value of "byCategory" separately, average
    byCategoryValues <- clusterSummaryTable[, unique(get(byCategory))]
    matrixSize <- clusterSummaryTable[, uniqueN(get(category))]
    corrMatrix <- matrix(0, nrow = matrixSize, ncol = matrixSize)
    for (byCategoryValue in byCategoryValues) {
      categoryValueTable <- dcast(clusterSummaryTable[get(byCategory) == byCategoryValue],
          value.var = "Value", formula = paste0("Setting ~ ", category))
      categoryValueTable[, Setting := NULL]
      # Ignore warning if standard deviation is zero
      corrMatrix <- corrMatrix + 1/length(byCategoryValues) *
        suppressWarnings(cor(categoryValueTable, use = "pairwise.complete.obs",
                             method = correlationMethod))
    }
  }
  if (plot) {
    if (is.na(plotFilePath)) {
      corrplot(corr = corrMatrix)
    } else {
      png(filename = paste0(plotFilePath, ".png"), width = 960, height = 540)
      corrplot(cor = corrMatrix)
      dev.off()
    }
  }
  return(corrMatrix)
}

# Goes through each directory in "dataPaths", finds the evaluation tables which have
# multiples values for the columns mentioned in "categories" (default: look at all
# experimental settings + indices + references), and saves the correlation plots
# between these experimental settings/indices in the corresponding plot directory
# (or the directories you provide for plotting).
# "infStrategy": see summarizeExperimentSettingCorrelation()
saveAllEvaluationCorrelationPlots <- function(dataPaths, plotPaths = NULL,
    categories = NULL, infStrategy = "drop") {
  for (i in seq_len(length(dataPaths))) {
    dataPath <- dataPaths[i]
    cat("Creating correlation plots for \"", dataPath, "\".\n", sep = "")
    if (is.null(plotPaths)) {
      plotPath <- paste0(gsub("data", "plots", dataPath), "ExperimentCorrelationPlots/")
    } else {
      plotPath <- plotPaths[i]
    }
    dir.create(plotPath, showWarnings = FALSE, recursive = TRUE)
    resultPath <- paste0(addSlashesToDirs(dataPath), "results/")
    # Only files which have multiple indices
    evaluationFilePaths <- list.files(path = resultPath, full.names = TRUE,
        pattern = "(clusterCount|entropy|internal|external|groundTruth|Forecast).*\\.rds$")
    for (evaluationFilePath in evaluationFilePaths) {
      cat("-File: \"", basename(evaluationFilePath), "\".-\n", sep = "")
      evaluationTable <- readRDS(evaluationFilePath)
      isForecastTable <- grepl("ForecastTable", evaluationFilePath)
      tableName <- gsub("\\.rds$", "", basename(evaluationFilePath))
      if (grepl("groundTruth", tableName)) {
        tableName <- paste0(strsplit(tableName, "_")[[1]][c(1,3)], collapse = "_")
      }
      if (is.null(categories)) { # take all experimental settings + Index + Reference
        relevantCols <- setdiff(colnames(evaluationTable), "Value")
      } else {
        relevantCols <- intersect(colnames(evaluationTable), categories)
      }
      for (colName in relevantCols) {
        if (evaluationTable[, uniqueN(get(colName))] == 1) {
          next # nothing to correlate
        }
        if (colName == "Reference" & !isForecastTable) {
          next  # continue; for forecasts, there is only one Index, but the different
          # staretegies are saved in the Reference column; for all other tables,
          # references should be analysed independently
        }
        if (!isForecastTable && "Reference" %in% names(evaluationTable)) {
          # Separate tables for references (affects internal and external CVIs)
          for (reference in evaluationTable[, unique(Reference)]) {
            summarizeExperimentSettingCorrelation(
              clusterSummaryTable = evaluationTable[Reference == reference],
              category = colName, infStrategy = infStrategy, plot = TRUE,
              plotFilePath = paste0(plotPath, "Correlation_", tableName,
                                    "_", reference, "_", colName))
          }
        } else {
          summarizeExperimentSettingCorrelation(clusterSummaryTable = evaluationTable,
              category = colName, infStrategy = infStrategy, plot = TRUE,
              plotFilePath = paste0(plotPath, "Correlation_", tableName,
                                    "_", colName))
        }
      }
    }
  }
  return(TRUE)
}

# Creates one plot summarizing an evaluation index with the function "funcName",
# grouping by "category"; saves plot under the "plotPath"
saveExperimentCategoryPlot <- function(clusterSummaryTable, funcName = "mean",
    plotName = "Evaluation", plotPath, category) {
  png(filename = paste0(plotPath, plotName, "_", funcName, "By", category, ".png"),
      width = 960, height = 540)
  plotSummaryTable <- clusterSummaryTable[,
      .(Value = do.call(funcName, list(Value, na.rm = TRUE))),
      by = c(category, "Level")]
  plotSummaryTable <- plotSummaryTable[!is.na(Value)]
  if (category == "Reference") {
    # Reference determines line type in ordinary evaluation, but we want it as
    # group (color) attribute
    setnames(plotSummaryTable, "Reference", "Strategy")
    category <- "Strategy"
  }
  print(plotOneClusteringEvaluation(plotSummaryTable, groupExpression = category))
  dev.off()
  return(TRUE)
}

# Creates plots summarizing an evaluation index with the function "funcName" by
# dissimilarity, aggregation and clustering (each combined with the aggregation
# levels); saves plots under the "plotPath"; you might want to limit your cluster
# summary table to one Index and one Reference if it has multiple ones
saveExperimentCategoryPlots <- function(clusterSummaryTable, funcName = "mean",
    plotName = "Evaluation", plotPath) {
  for (category in c("Dataset", "Dissimilarity", "Clustering", "Aggregation")) {
    if (clusterSummaryTable[, uniqueN(get(category))] > 1) {
      saveExperimentCategoryPlot(clusterSummaryTable = clusterSummaryTable,
          funcName = funcName, plotName = plotName, plotPath = plotPath,
          category = category)
    }
  }
  return(TRUE)
}

# Saves plots for an evaluation table comparing aggregation levels, considering 
# the values of the evaluation index ungrouped (boxplot) as well as grouped (using
# a summary statistic "funcName")
saveExperimentSummaryPlots <- function(clusterSummaryTable, funcName = "mean",
                                      plotName = "Overview", plotPath) {
  cat("Creating plots for ", plotName, ".\n", sep = "")
  dir.create(plotPath, showWarnings = FALSE, recursive = TRUE)
  # All values plotted against the aggregation level
  png(filename = paste0(plotPath, plotName, ".png"), width = 960, height = 540)
  print(plotClusteringEvaluationOverview(clusterSummaryTable, plotTitle = plotName))
  dev.off()
  # Comparison of categories for the aggregation levels
  saveExperimentCategoryPlots(clusterSummaryTable = clusterSummaryTable,
      funcName = funcName, plotName = plotName, plotPath = plotPath)
}

# Prints the standard R summary over all experiment categories
summarizeOverall <- function(clusterSummaryTable, hint = "") {
  cat("### Summary")
  if (hint != "") cat(" (", hint, ")", sep = "")
  cat("\n\n")
  print(clusterSummaryTable[, summary(Value)])
  cat("\n")
  return(TRUE)
}

# Determines the number of time series and machines from batch experiment result directories
summarizeTimeSeriesCount <- function(dataPaths, print = TRUE) {
  result <- data.table(Dataset = character(), DatasetSize = integer(), MachineCount = integer())
  for (dataPath in addSlashesToDirs(dataPaths)) {
    if (print) cat("Directory: ", basename(dataPath), "\n", sep = "")
    datasetName <- strsplit(list.files(paste0(dataPath, "aggData/"), pattern = "\\.rds")[1],
                            split = "_")[[1]][1]
    if (file.exists(paste0(dataPath, "results/machineSensorDayCombinations.rds"))) {
      # Fixed length, no forecasting
      overviewTable <- readRDS(paste0(dataPath, "results/machineSensorDayCombinations.rds"))
      datasetSize <- nrow(overviewTable)
      machineCount <- overviewTable[, uniqueN(machine_name)]
      if (print) cat("Time series: ", datasetSize, " - Machines: ", machineCount, "\n", sep = "")
      result <- rbind(result, list(Dataset = datasetName, DatasetSize = datasetSize,
                                   MachineCount = machineCount))
    } else if (file.exists(paste0(dataPath, "results/trainMachineSensorDayCombinations.rds"))) {
      # Fixed length, forecasting
      trainOverviewTable <- readRDS(paste0(dataPath, "results/trainMachineSensorDayCombinations.rds"))
      trainTotalLength <- nrow(trainOverviewTable)
      trainMachineCount <- trainOverviewTable[, uniqueN(machine_name)]
      testOverviewTable <- readRDS(paste0(dataPath, "results/testMachineSensorDayCombinations.rds"))
      testTotalLength <- nrow(testOverviewTable)
      testMachineCount <- testOverviewTable[, uniqueN(machine_name)]
      trainClusLength <- length(readRDS(paste0(dataPath,
          "results/trainDataPrevActList.rds"))[[1]]$prev)
      testClusLength <- length(readRDS(paste0(dataPath,
          "results/testDataPrevActList.rds"))[[1]]$prev)
      if (print) cat("Clustering time series: ", trainClusLength + testClusLength, " (",
          trainClusLength, ":", testClusLength, ") - Machines: ", trainMachineCount, ":",
          testMachineCount, "\n", sep = "")
      result <- rbind(result, list(Dataset = datasetName, DatasetSize = trainClusLength,
                                   MachineCount = trainMachineCount))
      if (print) cat("Dataset time series: ", trainTotalLength + testTotalLength, " (",
          trainTotalLength, ":", testTotalLength, ") - Machines: ", trainMachineCount, ":",
          testMachineCount, "\n", sep = "")
    } else if (file.exists(paste0(dataPath, "results/timeSeriesOverviewTable.rds"))) {
      # Variable length
      overviewTable <- readRDS(paste0(dataPath, "results/timeSeriesOverviewTable.rds"))
      datasetSize <- nrow(overviewTable)
      machineCount <- overviewTable[, uniqueN(machine_name)]
      if (print) cat("Time series:", datasetSize, " - Machines: ", machineCount, "\n", sep = "")
      result <- rbind(result, list(Dataset = datasetName, DatasetSize = datasetSize,
                                   MachineCount = machineCount))
    } else {
      stop("Expected result files not found at the specified path.")
    }
    if (print) cat("\n")
  }
  return(result)
}

# Adds slashes to the end of paths if not existing
addSlashesToDirs <- function(dirs) {
  return(sapply(dirs, function(dirName) {
    if (substr(dirName, start = nchar(dirName), stop = nchar(dirName)) != "/") {
      return(paste0(dirName, "/"))
    } else {
      return(dirName)
    }
  }, USE.NAMES = FALSE))
}

# Merges evaluation tables from multiple datasets
createMergedEvaluationTable <- function(dataPaths, fileName) {
  result <- rbindlist(lapply(addSlashesToDirs(dataPaths), function(dataPath) {
    return(readRDS(paste0(dataPath, "results/", fileName)))
  }))
  result <- formatAndSortClusterSummaryTable(result)
  return(result)
}

# Merges all evaluation tables from the directories in the vector "dataPaths" and
# writes the results to "ouputPath".
saveAllMergedEvaluationTables <- function(dataPaths, outputPath) {
  outputPath <- addSlashesToDirs(outputPath)
  if (!grepl("results/$", outputPath)) {
    outputPath <- paste0(outputPath, "results/")
  }
  dir.create(outputPath, showWarnings = FALSE, recursive = TRUE)
  resultPaths <- paste0(addSlashesToDirs(dataPaths), "results/")
  # Not all files are suitable for merge, e.g. params or time series defintions
  evaluationFilesPathList <- lapply(resultPaths, list.files, full.names = TRUE,
      pattern = "(Time|clusterCount|entropy|internal|external|groundTruth|Forecast).*\\.rds$")
  evaluationFilesList <- lapply(evaluationFilesPathList, basename)
  # Ground truth tables contain dataset names which we remove before checking
  # if same files exist for all datasets
  datasetRenamedFilesLists <- lapply(evaluationFilesList, function(evaluationFiles) {
    sapply(evaluationFiles, function(fileName) {
      if (grepl("_", fileName)) {
        newFileName <- rep("DatasetsMerged", 3)
        newFileName[c(1,3)] <- strsplit(fileName, "_")[[1]][c(1,3)]
        return(paste0(newFileName, collapse = "_"))
      } else {
        return(fileName)
      }
    }, USE.NAMES = FALSE)
  })
  for (i in 2:length(resultPaths)) {
    if (length(datasetRenamedFilesLists[[1]]) != length(datasetRenamedFilesLists[[i]]) ||
        any(datasetRenamedFilesLists[[1]] != datasetRenamedFilesLists[[i]])) {
      cat("[WARN] Directory \"", resultPaths[1], "\" and directory \"", resultPaths[i],
          "\" contain different evaluation files.\n", sep = "")
    }
  }
  for (evaluationFileName in unique(unlist(datasetRenamedFilesLists))) {
    cat("Merging \"", evaluationFileName, "\" files.\n", sep = "")
    mergedData <- rbindlist(lapply(evaluationFilesPathList, function(evalFilePaths) {
      if (grepl("groundTruth", evaluationFileName)) {
        # remove dataset name and check for 1st + last part of name
        fileNameParts <- strsplit(evaluationFileName, "_")[[1]][c(1,3)]
        matchingFile <- intersect(grep(pattern = paste0(fileNameParts[1], "_"),
                x = evalFilePaths, value = TRUE, fixed = TRUE),
            grep(pattern = paste0("_", fileNameParts[2]), x = evalFilePaths,
                 value = TRUE, fixed = TRUE))
      } else {
        matchingFile <- grep(pattern = evaluationFileName, x = evalFilePaths,
                             value = TRUE, fixed = TRUE)
      }
      if (length(matchingFile) == 0) {
        return(NULL)
      } else if (length(matchingFile) == 1) {
        return(readRDS(matchingFile))
      } else {
        stop("Error in file matching procedure.")
      }
    }))
    mergedData <- formatAndSortClusterSummaryTable(mergedData)
    saveRDS(mergedData, file = paste0(outputPath, evaluationFileName))
  }
  return(TRUE)
}

# Merged the cluster count evaluation table with other tables and renames
# columns; assumes a cluster count and an entropy table merged over several
# datasets have already been created ("mergedDataPath"), but acesses the
# original, one-dataset data as well ("inputDataPaths") to determine the size
# of the biggest clusters
prepareClusterCountTableForPlots <- function(mergedDataPath, inputDataPaths) {
  clusterCountTable <- readRDS(paste0(mergedDataPath, "results/clusterCountTable.rds"))
  setnames(clusterCountTable, "Value", "ClusterCount")
  # Consider dataset size
  timeSeriesCounts <- summarizeTimeSeriesCount(inputDataPaths, print = FALSE)
  clusterCountTable <- merge(clusterCountTable, timeSeriesCounts)
  clusterCountTable[, RelClusterCount := ClusterCount / DatasetSize]
  entropyTable <- readRDS(paste0(mergedDataPath, "results/entropyTable.rds"))
  # Merge with entropy table
  setnames(entropyTable, "Value", "ClusterSizeEntropy")
  clusterCountTable <- formatAndSortClusterSummaryTable(clusterCountTable)
  clusterCountTable <- merge(clusterCountTable, entropyTable)
  # Consider size of biggest cluster (not part of standard evaluation table)
  maxClusterSizeTable <- rbindlist(lapply(inputDataPaths, function(dataDir) {
    createMaxClusterSizeTable(readRDS(paste0(dataDir, "results/clusteringResult.rds")))
  }))
  setnames(maxClusterSizeTable, "Value", "MaxClusterSize")
  clusterCountTable <- merge(clusterCountTable, maxClusterSizeTable)
  clusterCountTable[, RelMaxClusterSize := MaxClusterSize / DatasetSize]
  return(clusterCountTable)
}

# Merges two tables, one containing experiments with fixed length time series,
# the other one containing experiments with variable length time series
# Only the dissimilarities present in both tables are used.
mergeFixedAndVariableLengthTables <- function(fLengthSummaryTable, vLengthSummaryTable) {
  commonDissimilarities <- intersect(fLengthSummaryTable[, unique(Dissimilarity)],
                                     vLengthSummaryTable[, unique(Dissimilarity)])
  fLengthSummaryTable <- copy(fLengthSummaryTable)
  fLengthSummaryTable <- fLengthSummaryTable[Dissimilarity %in% commonDissimilarities]
  fLengthSummaryTable[, Dataset := paste(Dataset, "(F)")]
  vLengthSummaryTable <- copy(vLengthSummaryTable)
  vLengthSummaryTable <- vLengthSummaryTable[Dissimilarity %in% commonDissimilarities]
  vLengthSummaryTable[, Dataset := paste(Dataset, "(V)")]
  if (nrow(fLengthSummaryTable) != nrow(vLengthSummaryTable)) {
    warning(paste0("The two tables seem to contain different experimental settings ",
                   "apart from the dissimilarities."))
  }
  result <- rbind(fLengthSummaryTable, vLengthSummaryTable)
  result <- formatAndSortClusterSummaryTable(result)
  return(result)
}
