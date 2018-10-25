# Call with "Rscript batchscripts/BatchExperimentReportScript.R <<dataPath>>"
# from project directory

source("BatchExperimentEvaluationUtility.R")

### Global params ###

dataPath <- commandArgs(trailingOnly = TRUE)[1]

### Read experiment data ###

cat("---Loading data---\n")
resultsPath <- paste0(dataPath, "results/")
# Read all *.rds file in the "results" directory
for (resultsFileName in list.files(resultsPath, full.names = TRUE, pattern = "\\.rds")) {
  assign(gsub(pattern = ".rds", replacement = "", x = basename(resultsFileName), fixed = TRUE),
         value = readRDS(resultsFileName))
}

### Evaluation ###

cat("---Writing evaluation to file---\n")
sink(file = paste0(resultsPath, "evaluation.txt"))
cat("# Evaluation for \"", dataPath, "\" (created at ", format(Sys.time()), ")\n\n", sep = "")

## Cluster plots ##
# Example, you need to decide on one clustering/aggregation/level/dissimilarity
# saveExperimentClusterPlots(dataPath = dataPath, clusterAssignmentsTable = clusteringResult,
#     aggregation = "mean", level = 1440, dissimilarity = "L2", clustering = "PAM")

## Cluster counts ##
cat("## Cluster counts\n\n")
invisible(summarizeOverall(clusterCountTable))
invisible(summarizeExperimentCategories(clusterCountTable, funcName = "mean"))

## Entropy ##
cat("## Entropy\n\n")
invisible(summarizeOverall(entropyTable))
invisible(summarizeExperimentCategories(entropyTable, funcName = "mean"))
cat("### NA check\n\n")
invisible(checkEntropyNAs(entropyTable))
cat("\n")

## Internal Validity ##
index <- "Silhouette"
cat("## Internal CVIs (", index, ")\n\n", sep = "")
invisible(summarizeOverall(internalCVITable[Index == index & Reference == "B"], hint = "base"))
invisible(summarizeOverall(internalCVITable[Index == index & Reference == "C"], hint = "current"))
invisible(summarizeExperimentCategories(internalCVITable[Index == index & Reference == "B"],
                                        funcName = "mean", hint = "base"))
invisible(summarizeExperimentCategories(internalCVITable[Index == index & Reference == "C"],
                                        funcName = "mean", hint = "current"))
if (exists("clusteringResult")) { # might not exist in a multi-dataset setting
  cat("### NA check\n\n")
  invisible(checkInternalCVINAs(internalCVITable, clusterAssignmentsTable = clusteringResult))
  cat("\n")
}

## External validity ##
index <- grep("Dongen", externalCVITable[, unique(Index)], value = TRUE)
cat("## External CVIs (", index, ")\n\n", sep = "")
invisible(summarizeOverall(externalCVITable[Index == index & Reference == "B"], hint = "base"))
invisible(summarizeOverall(externalCVITable[Index == index & Reference == "P"], hint = "previous"))
invisible(summarizeExperimentCategories(externalCVITable[Index == index & Reference == "B"],
                                        funcName = "mean", hint = "base"))
invisible(summarizeExperimentCategories(externalCVITable[Index == index & Reference == "P"],
                                        funcName = "mean", hint = "previous"))
cat("### NA check\n\n")
invisible(checkExternalCVINAs(externalCVITable, clusterCountTable = clusterCountTable))
cat("\n")

## Ground truth ##
for (gtTableName in ls(pattern = "groundTruth")) {
  cat("## Ground truth (", gtTableName, ")\n\n", sep = "")
  currentGroundTruthTable <- get(gtTableName)
  # Be aware that purity/uniformity of clusters tends to be high if there are
  # many/small clusters (and a few classes or even one dominating class) and
  # purity/uniformity of classes tends to be high if there are many/small classes
  # and/or a few clusters; van Dongen/NMI balance this)
  cat("### Number of values > .7 per index\n")
  print(currentGroundTruthTable[Value > .7, table(Index)])
  cat("\n")
  cat("### Correlation between cluster count and index\n\n")
  # warning if sd() is 0, but returns NA anyway -> suppress warning
  print(currentGroundTruthTable[, suppressWarnings(cor(Value, clusterCountTable$Value,
                                                       use = "na.or.complete")), by = Index])
  cat("\n")
  index <- grep("Dongen", externalCVITable[, unique(Index)], value = TRUE)
  invisible(summarizeOverall(currentGroundTruthTable[Index == index], hint = index))
  invisible(summarizeExperimentCategories(currentGroundTruthTable[Index == index],
                                          funcName = "mean", hint = index))
  invisible(tableExperimentCategories(currentGroundTruthTable[Index == index & Value > .7],
                                      hint = paste0(index, ", Value > 0.7")))
  cat("### NA check\n\n")
  invisible(checkGroundTruthNAs(currentGroundTruthTable))
  cat("\n")
  rm(currentGroundTruthTable)
}

## Forecasting ## (not in all experiments)
for (forecastTableName in ls(pattern = "ForecastTable")) {
  cat("## Forecasting error ratio (", forecastTableName, ")\n\n", sep = "")
  forecastErrorTable <- get(forecastTableName)
  forecastErrorRatioTable <- transformToErrorRatioTable(forecastErrorTable)
  cat("### Mean per technique\n\n")
  print(forecastErrorRatioTable[, mean(Value), by = Reference])
  cat("\n")
  invisible(summarizeExperimentCategories(forecastErrorRatioTable[Reference == "clus.weighted"],
                                          funcName = "mean", hint = "clustering-based (weighted) forecast"))
  invisible(tableExperimentCategories(forecastErrorRatioTable[Reference == "clus.weighted" & Value < .9],
                                      hint = "clustering-based (weighted) forecast, Value < 0.9"))
  cat("### NA check\n\n")
  invisible(checkForecastErrorNAs(forecastErrorTable))
  cat("\n")
  rm(forecastErrorTable)
  rm(forecastErrorRatioTable)
}

## Forecast plots ##
# Example, you need to decide on one clustering/aggregation/level/dissimilarity
# saveExperimentForecastPlots(dataPath = dataPath, aggregation = "mean", level = 1440,
#     dissimilarity = "L2", clustering = "PAM")

sink()
cat("---Finished report generation---\n")
