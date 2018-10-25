# Call with "Rscript batchscripts/BatchExperimentPlotScript.R <<dataPath>>"
# from project directory

source("BatchExperimentEvaluationUtility.R")

### Global params ###

dataPath <- commandArgs(trailingOnly = TRUE)[1]
plotPath <- paste0(gsub("data", "plots", dataPath), "ExperimentSummaryPlots/")
correlationPlotPath <- paste0(gsub("data", "plots", dataPath), "ExperimentCorrelationPlots/")
AGG_FUNC <- "mean"

### Read experiment data ###

cat("---Loading data---\n")
resultsPath <- paste0(dataPath, "results/")
# Read all *.rds file in the "results" directory
for (resultsFileName in list.files(resultsPath, full.names = TRUE, pattern = "\\.rds")) {
  assign(gsub(pattern = ".rds", replacement = "", x = basename(resultsFileName), fixed = TRUE),
         value = readRDS(resultsFileName))
}

### Evaluation ###

cat("---Creating evaluation plots---\n")

invisible(saveExperimentSummaryPlots(clusterCountTable, funcName = AGG_FUNC,
    plotName = "ClusterCount", plotPath = plotPath))
invisible(saveExperimentSummaryPlots(entropyTable, funcName = AGG_FUNC,
    plotName = "ClusterSizeEntropy", plotPath = plotPath))
invisible(saveExperimentSummaryPlots(internalCVITable[Index == "Silhouette" &
    Reference == "B"], funcName = AGG_FUNC, "Silhouette_Base", plotPath = plotPath))
invisible(saveExperimentSummaryPlots(internalCVITable[Index == "Silhouette" &
    Reference == "C"], funcName = AGG_FUNC, "Silhouette_Current", plotPath = plotPath))
defaultExternalIndex <- grep("Dongen", externalCVITable[, unique(Index)], value = TRUE)
invisible(saveExperimentSummaryPlots(externalCVITable[Index == defaultExternalIndex &
    Reference == "B"], funcName = AGG_FUNC, paste0(defaultExternalIndex, "_Base"),
    plotPath = plotPath))
invisible(saveExperimentSummaryPlots(externalCVITable[Index == defaultExternalIndex &
    Reference == "P"], funcName = AGG_FUNC, paste0(defaultExternalIndex, "_Previous"),
    plotPath = plotPath))
for (gtTableName in ls(pattern = "groundTruth")) {
  currentGroundTruthTable <- get(gtTableName)
  gtTableName <- strsplit(gtTableName, "_")[[1]][3]
  invisible(saveExperimentSummaryPlots(currentGroundTruthTable[Index == defaultExternalIndex],
      funcName = AGG_FUNC, paste0(defaultExternalIndex, "_", gtTableName), plotPath = plotPath))
  rm(currentGroundTruthTable)
}
for (forecastTableName in ls(pattern = "ForecastTable")) {
  forecastErrorTable <- get(forecastTableName)
  forecastTableName <- gsub("Table", "", forecastTableName)
  forecastErrorRatioTable <- transformToErrorRatioTable(forecastErrorTable)
  invisible(saveExperimentSummaryPlots(forecastErrorRatioTable[Reference == "clus.weighted"],
      funcName = AGG_FUNC,plotName = paste0(forecastTableName, "_clusWeighted"),
      plotPath = plotPath))
  saveExperimentCategoryPlot(forecastErrorRatioTable, funcName = AGG_FUNC,
      plotName = forecastTableName, plotPath = plotPath, category = "Reference")
  saveExperimentCategoryPlot(forecastErrorRatioTable[Dissimilarity == "L2"],
      funcName = AGG_FUNC, plotName = paste0(forecastTableName, "_L2"), plotPath = plotPath,
      category = "Reference")
}

# Correlation plot routine automatically goes through all tables
saveAllEvaluationCorrelationPlots(dataPaths = dataPath,
    plotPaths = correlationPlotPath, infStrategy = "drop")

cat("---Finished plotting---\n")
