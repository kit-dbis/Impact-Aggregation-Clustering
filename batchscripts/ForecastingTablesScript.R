# Simplified script if you just want to (re-)generate the forecast tables (with existing
# clustering and assignments); use: "Rscript batchscripts/ForecastingTablesScript.R <<dataPath>"
source("BatchUtility.R")

startTime <- Sys.time()
cat("-----Forecasting table routine started at ", format(startTime), "-----\n", sep = "")

# Parameter initialization
dataPath <- commandArgs(trailingOnly = TRUE)[1]
resultsPath <- paste0(dataPath, "results/")
sink(file = paste0(resultsPath, "params.txt"), append = TRUE)
cat("\nForecasting tables re-created at ", format(startTime), ", commit: ",
    system("git log --oneline -1", intern = TRUE), "\n", sep = "")
sink()
datasetName <- readRDS(paste0(resultsPath, "params.rds"))$datasetName
fullAggDataBasePathList <- list(paste0(dataPath, "fullTrainAggData/"))
names(fullAggDataBasePathList) <- datasetName
fullTestDataBasePathList <- list(paste0(dataPath, "fullTestAggData/"))
names(fullTestDataBasePathList) <- datasetName
trainDataPrevActList <- readRDS(paste0(resultsPath, "trainDataPrevActList.rds"))
trainDataAssignmentsTable <- readRDS(paste0(resultsPath, "clusteringResult.rds"))
testDataPrevActList <- readRDS(paste0(resultsPath, "testDataPrevActList.rds"))
testDataDissMatrixPath <- paste0(dataPath, "testDissMatrices/")

# Forecasting: copied + adapted from BatchUtility::saveForecastEvaluationTables()
cat("---Creating forecast evaluation tables---\n")
# Train data input + target series
fullTrainDataList <- list()
for (datasetName in names(fullAggDataBasePathList)) {
  fullTrainDataList[[datasetName]] <- readRDS(list.files(
    fullAggDataBasePathList[[datasetName]], pattern = "_Raw_.*\\.rds",
    full.names = TRUE)[1])
}
# Test data input + target series
fullTestDataList <- list()
for (datasetName in names(fullTestDataBasePathList)) {
  fullTestDataList[[datasetName]] <- readRDS(list.files(
    fullTestDataBasePathList[[datasetName]], pattern = "_Raw_.*\\.rds",
    full.names = TRUE)[1])
}
# Load assignment from disk
testDataAssignmentResult <- readRDS(paste0(resultsPath, "testDataAssignmentResult.rds"))
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
saveRDS(baseForecastTable, file = paste0(resultsPath, "baseForecastTable.rds"))
# Forecast - aggregated data (data used for the clustering on each level)
cat("--Creating relative forecast evaluation table--\n")
startTime <- as.numeric(Sys.time())
relativeForecastTable <- createRelativeForecastErrorTable(
  newDatasetBasePathList = fullTestDataBasePathList,
  oldDatasetBasePathList = fullAggDataBasePathList,
  newDataPrevActList = testDataPrevActList, oldDataPrevActList = trainDataPrevActList,
  newDataAssignmentsTable = testDataAssignmentResult,
  oldDataAssignmentsTable = trainDataAssignmentsTable,
  dissMatrixPaths = testDataDissMatrixPath, dissWeightConversion = "gaussian",
  forecastAggregation = "mean"
)
endTime <- as.numeric(Sys.time())
cat("Creating relative forecast evaluation table took ", round(endTime - startTime), "s\n", sep = "")
saveRDS(relativeForecastTable, file = paste0(resultsPath, "relativeForecastTable.rds"))

cat("-----Forecasting table routine finished at ", format(Sys.time()), "-----\n", sep = "")
