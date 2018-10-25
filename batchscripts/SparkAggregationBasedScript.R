# Call this script with "Rscript batchscripts/SparkAggregationBasedScript.R <<arguments>>"
# from project dir
# Supply the following arguments:
# --dataPath="path/to/dir/" : directory of input data (will also be used to save
# output); input are CSV files with aggregations exported from Spark
# --plotPath="path/to/dir/" : directory to place cluster plots (optional)
# --datasetName="YourDatasetName" : will be used for time series lists and in
# evaluation tables; !!!must not contain an underscore!!!
# --normalize : if time series should be z-scored (individually for each
# aggregation, level, series) when creating them (affects all further steps of
# analysis); leave parameter out if you want no normalization
# --testRatio=0.5 : (relative) number (in [0,1)) of time series used as
# test set in forecasting (optional, default is 0)
# --aggregations="mean|min_max|median" : aggregations to be used, multivariate
# aggregations denoted by underscore (_), aggregation types by pipe (|)
# --machines="AVT_01|AVT_02" : machines to be used (default:all), multiple
# machines separated by pipe character
# --intervalLength=86400 : maximum length (in seconds) of the time interval
# represented by one time series (default: one day, but you might also want to
# choose a divisor, e.g. 3600 for time series representing hours); can also be
# negative if you want variable length time series; can be zero if you want
# variable length time series and the length is determined by a separate *.rds
# file containing machine, start and end time

source("BatchUtility.R")

cat("-----Experiment started at ", format(Sys.time()), "-----\n", sep = "")

# Input parameter handling
params <- parseInputParams()

# Pre-processing (if not already pre-processed from previous run)
aggDataPath <- paste0(params$dataPath, "aggData/")
resultsPath <- paste0(params$dataPath, "results/")
testAggDataPath <- paste0(params$dataPath, "testAggData/")
fullForecastAggDataPath <- paste0(params$dataPath, "fullTrainAggData/")
fullForecastTestAggDataPath <- paste0(params$dataPath, "fullTestAggData/")
aggNamesList <- convertAggregationStringToList(params$aggregations)
aggFunctionsList <- sparkAggNameListToAggFunctionList(aggNamesList,
    aggregationLevels = getDefaultDailyAggregationLevels(params$intervalLength),
    variableLength = params$intervalLength <= 0)
savePreprocessedTimeSeries(
  datasetName = params$datasetName,
  dataPath = params$dataPath,
  aggDataPath = aggDataPath,
  resultsPath = resultsPath,
  machineNames = convertMachineStringToVector(params$machines),
  intervalLength = params$intervalLength,
  aggNamesList = aggNamesList,
  aggFunctionsList = aggFunctionsList,
  normalize = params$normalize,
  testAggDataPath = testAggDataPath, # forecasting is optional, might not be used
  fullForecastAggDataPath = fullForecastAggDataPath,
  fullForecastTestAggDataPath = fullForecastTestAggDataPath,
  testRatio = params$testRatio
)

# Clustering, evaluation - entropy, internal, external, ground truth
aggDataBasePathList <- list(aggDataPath)
names(aggDataBasePathList) <- params$datasetName
dissMatrixPath <- paste0(params$dataPath, "dissMatrices/")
dissMatrixPathList <- list(dissMatrixPath)
names(dissMatrixPathList) <- params$datasetName
if (params$plotPath == "") {
  clusPlotBasePathList <- list()
} else {
  clusPlotBasePathList <- list(params$plotPath)
  names(clusPlotBasePathList) <- params$datasetName
}
disFunctionsList <- getDefaultDissimilarities(aggNamesList, params$intervalLength)
runImpactAnalysis(
  datasetList = NULL, # already aggregated
  aggDataBasePathList = aggDataBasePathList,
  dissMatrixPathList = dissMatrixPathList,
  tsPlotPathList = list(), # no time series plots to copy
  clusPlotBasePathList = clusPlotBasePathList,
  resultsPath = resultsPath,
  aggList = aggFunctionsList,
  disList = disFunctionsList,
  clusList = DEFAULT_CLUSTERINGS
)

# Evaluation - forecasting (optional)
if (params$testRatio > 0) {
  fullAggDataBasePathList <- list(fullForecastAggDataPath)
  names(fullAggDataBasePathList) <- params$datasetName
  testDataBasePathList <- list(testAggDataPath)
  names(testDataBasePathList) <- params$datasetName
  fullTestDataBasePathList <- list(fullForecastTestAggDataPath)
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
    testDataDissMatrixPath = paste0(params$dataPath, "testDissMatrices/"),
    aggList = aggFunctionsList, disList = disFunctionsList, outputPath = resultsPath
  )
}

cat("-----Experiment finished at ", format(Sys.time()), "-----\n", sep = "")

cat("-----Running report script------\n")
system(paste0("Rscript batchscripts/BatchExperimentReportScript.R ", params$dataPath))

cat("-----Running plot script------\n")
system(paste0("Rscript batchscripts/BatchExperimentPlotScript.R ", params$dataPath))
