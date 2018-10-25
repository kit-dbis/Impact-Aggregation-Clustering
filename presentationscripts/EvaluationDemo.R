source("BatchExperimentEvaluationUtility.R")

### Global parameters ###

# With "experiment", we refer to one execution of the main experimentation script.
# If you execute that script multiple times, i.e., for different datasets, each
# experiment should go into a different sub-directory. Some of the routines here
# combine the results from multiple experiments.
DATA_DIR <- "<<dir/with/several/experiments/>"
forecastDirs <- "<<dir/with/several/experiments/>" # should contain (several) experiment dirs with forecasting evaluation enabled
experimentDirs <- "<<dir/with/several/experiments/>" # should contain (several) experiment dirs

############################################################

### Check number of time series and plot time series for one dataset ###

dataPath = "<<dir/of/one/experiment/>>"
summarizeTimeSeriesCount(dataPath)
saveExperimentBaseTimeSeriesPlots(dataPath)
# Plot selected clusterings
saveExperimentClusterPlots(dataPath = dataPath, aggregation = "mean",
    level = c(24, 144, 1440), dissimilarity = c("L2", "DTW"),
    clustering = c("PAM", "Hier.avg", "DBSCAN", "AffProp"))
# If forecasting used: plot selected forecasting results
saveExperimentForecastPlots(dataPath = dataPath, aggregation = "mean",
    level = c(24, 144, 1440), dissimilarity = c("L2", "DTW"),
    clustering = c("PAM", "Hier.avg", "DBSCAN", "AffProp"),
    relative = c(FALSE, TRUE), withBackground = FALSE)

# Evaluate aggregation trade-off (example: validity vs execution time)

# Obtain input tables from any dataset (or merged over datasets)
externalCVITable <- readRDS(paste0(dataPath, "results/externalCVITable.rds"))
internalCVITable <- readRDS(paste0(dataPath, "results/internalCVITable.rds"))
clusteringTimeTable <- readRDS(paste0(dataPath, "results/clusteringTimeTable.rds"))
# Prepare external CVI table: bring back the "Raw" aggregation which is removed
# in the creation of the table to make more simple Level comparisons
externalCVITable <- externalCVITable[Reference == "B" & Index == "inv.vanDongen"]
externalCVITable <- revertRawResultsInClusterSummaryTable(externalCVITable)
externalCVITable[, c("Index", "Reference") := NULL]
setnames(externalCVITable, "Value", "ExternalQuality")
# Prepare internal CVI table: one index, map to [0,1]
internalCVITable <- internalCVITable[Reference == "B" & Index == "Silhouette"]
internalCVITable[, Value := (Value + 1) / 2]
internalCVITable[, c("Index", "Reference") := NULL]
setnames(internalCVITable, "Value", "InternalQuality")
# Prepare computation time table
clusteringTimeTable[, TotalTime := ClusteringTime + DissimilarityTime]
clusteringTimeTable[, c("DissimilarityTime", "ClusteringTime") := NULL]
# Merge
tradeoffTable <- merge(externalCVITable, internalCVITable)
tradeoffTable <- merge(tradeoffTable, clusteringTimeTable)
# Compute Quality-per-Cost, format table for evaluation routines
tradeoffTable[, QpC := 0.5 * (InternalQuality + ExternalQuality) / TotalTime]
formattedTradeoffTable <- tradeoffTable[, .(Dataset, Dissimilarity, Clustering,
                                            Aggregation, Level, Value = QpC)]
# Find best experimental setting per category
summarizeBestExperimentSettings(clusterSummaryTable = formattedTradeoffTable,
    maximize = TRUE)
countBestExperimentSetting(formattedTradeoffTable[Dissimilarity == "DTW"],
    categoryName = "Level", relativeCount = TRUE)

### Compare any (optimizable) evaluation table over datasets ###

evaluationTable <- createMergedEvaluationTable(dataPaths = experimentDirs,
    fileName = "internalCVITable.rds")[Index != "Gen.Davies.Bouldin"]
# Ranking of experimental settings overall
summarizeBestExperimentSettings(evaluationTable, maximize = TRUE)
# Best experimental setting per dataset
summarizeBestExperimentSettings(evaluationTable, group = "Dataset", maximize = TRUE)
# Correlation analysis between indices
summarizeExperimentSettingCorrelation(evaluationTable, category = "Index")
# Save correlation plots for all tables (you can pass dimensions to compare as param)
saveAllEvaluationCorrelationPlots(dataPaths = experimentDirs, infStrategy = "drop")
# Plot comparing datasets at levels
plotOneClusteringEvaluation(evaluationTable[Index == "Silhouette" & Reference == "C",
    .(Value = mean(Value, na.rm = TRUE)), by = .(Dataset, Level)],
    groupExpression = "Dataset")
# Level boxplot over all datasets
plotClusteringEvaluationOverview(evaluationTable[Index == "Silhouette" & Reference == "B"])

### Compare evaluation measures ###

# make sure to limit to one Index, one Reference
internalCVITable <- createMergedEvaluationTable(dataPaths = experimentDirs,
    fileName = "internalCVITable.rds")
externalCVITable <- createMergedEvaluationTable(dataPaths = experimentDirs,
    fileName = "externalCVITable.rds")
evaluationTable1 <- internalCVITable[Reference == "B" & Level != 2880]
evaluationTable2 <- externalCVITable[Reference == "B" & Level != 2880]
evaluationTable1 <- formatClusterSummaryTableForJoin(evaluationTable1)
evaluationTable2 <- formatClusterSummaryTableForJoin(evaluationTable2)
checkDesignColumnIdentity(evaluationTable1, evaluationTable2)
# Single correlation
cor(evaluationTable1[Index == "Silhouette", Value],
    evaluationTable2[Index == "inv.vanDongen", Value], use = "na.or.complete")
# All correlations
evaluationTable <- rbind(evaluationTable1, evaluationTable2)
summarizeExperimentSettingCorrelation(evaluationTable, category = "Index",
                                      infStrategy = "drop")

### Compare forecasting over datasets ###

baseForecastErrorRatioTable <- transformToErrorRatioTable(
  createMergedEvaluationTable(dataPaths = forecastDirs, fileName = "baseForecastTable.rds"),
  baseRef = "naive")
relativeForecastErrorRatioTable <- transformToErrorRatioTable(
  createMergedEvaluationTable(dataPaths = forecastDirs, fileName = "relativeForecastTable.rds"),
  baseRef = "naive")
errorRatioTable <- baseForecastErrorRatioTable # choose one of the tables
# Which forecasting technique wins for each experimental settings? (lowest MAE)
countBestExperimentSetting(errorRatioTable, categoryName = "Reference", maximize = FALSE)
# Which forecasting technique looses for each experimental settings?
countBestExperimentSetting(errorRatioTable, categoryName = "Reference", maximize = TRUE)
# Which dataset has highest amount of wins for clustering-based forecast?
errorRatioTable[, countBestExperimentSetting(.SD, categoryName = "Reference", maximize = FALSE),
                by = Dataset][, N[Reference == "clus.weighted"] / sum(N), by = Dataset][order(-V1)]
# Which dataset has the lowest mean forecasting error for clustering-based forecast?
errorRatioTable[Reference == "clus.weighted", mean(Value), by = Dataset][order(V1)]

### Create merged evaluation tables ###

# Do not only merge one specific kind of tables, but all ones
saveAllMergedEvaluationTables(dataPaths = experimentDirs,
    outputPath = paste0(DATA_DIR, "Merged/"))
