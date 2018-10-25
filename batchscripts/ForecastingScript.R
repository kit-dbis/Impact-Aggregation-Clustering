# Shortened version of the main script which only performs test data assignment +
# forecasting (assuming pre-processing/clustering is completed and results are saved);
# usage: "Rscript batchscripts/ForecastingScript.R <<dataPath>"
source("BatchUtility.R")

startTime <- Sys.time()
cat("-----Forecasting routine started at ", format(startTime), "-----\n", sep = "")

# Parameter initialization
dataPath <- commandArgs(trailingOnly = TRUE)[1]
resultsPath <- paste0(dataPath, "results/")
sink(file = paste0(resultsPath, "params.txt"), append = TRUE)
cat("\nForecasting re-run at ", format(startTime), ", commit: ",
    system("git log --oneline -1", intern = TRUE), "\n", sep = "")
sink()
params <- readRDS(paste0(resultsPath, "params.rds"))
aggNamesList <- convertAggregationStringToList(params$aggregations)
aggFunctionsList <- sparkAggNameListToAggFunctionList(aggNamesList,
    aggregationLevels = getDefaultDailyAggregationLevels(params$intervalLength),
    variableLength = params$intervalLength <= 0)
disFunctionsList <- getDefaultDissimilarities(aggNamesList, params$intervalLength)
aggDataBasePathList <- list(paste0(dataPath, "aggData/"))
names(aggDataBasePathList) <- params$datasetName
fullAggDataBasePathList <- list(paste0(dataPath, "fullTrainAggData/"))
names(fullAggDataBasePathList) <- params$datasetName
testDataBasePathList <- list(paste0(dataPath, "testAggData/"))
names(testDataBasePathList) <- params$datasetName
fullTestDataBasePathList <- list(paste0(dataPath, "fullTestAggData/"))
names(fullTestDataBasePathList) <- params$datasetName

saveForecastEvaluationTables(
  fullTrainDataList = NULL, # already aggregated
  inputTrainDataBasePathList = aggDataBasePathList,
  fullTrainDataBasePathList = fullAggDataBasePathList,
  trainDataPrevActList = readRDS(paste0(resultsPath, "trainDataPrevActList.rds")), 
  trainDataAssignmentsTable = readRDS(paste0(resultsPath, "clusteringResult.rds")),
  fullTestDataList = NULL, # already aggregated
  inputTestDataBasePathList = testDataBasePathList,
  fullTestDataBasePathList = fullTestDataBasePathList,
  testDataPrevActList = readRDS(paste0(resultsPath, "testDataPrevActList.rds")),
  testDataDissMatrixPath = paste0(dataPath, "testDissMatrices/"),
  aggList = aggFunctionsList, disList = disFunctionsList, outputPath = resultsPath
)

cat("-----Forecasting routine finished at ", format(Sys.time()), "-----\n", sep = "")
