library(clv)
library(data.table)
library(doParallel)
library(FastTSDistances)
library(foreach)
library(ggplot2)
library(tcltk)

source("ClusteringUtility.R") # validity indices, checkAssignment

### Table transformations ###

# Returns a data.table containing experiment parameters (without assignments)
# including execution times; removes all "Time" columns from the table (in-place)
extractExecutionTimeTable <- function(clusterSummaryTable) {
  # instead of "noAssignmentColumns", we could also hard-code the columns we
  # want to keep in the executionTimeTable
  nonAssignmentColumns <- setdiff(colnames(clusterSummaryTable),
      grep(pattern = "^A[0-9]+", x = colnames(clusterSummaryTable), value = TRUE))
  timeColumns <- grep(pattern = "Time", x = colnames(clusterSummaryTable), value = TRUE)
  executionTimeTable <- clusterSummaryTable[, mget(nonAssignmentColumns)]
  clusterSummaryTable[, (timeColumns) := NULL]
  return(executionTimeTable)
}

# Takes a table with clustering results (assignments or any of the evaluation
# tables created with "create...Table()"), containing one Aggregation="Raw" row
# for each Clustering, Dissimilarity, Dataset and replaces them with multiple
# rows, one for each aggregation type (so we have base values for each aggregation
# and other evaluation methods can work with it).
# Furthermore, the column types are adapted and the table is sorted.
addRawResultsToClusterSummaryTable <- function(clusterSummaryTable) {
  aggRawTable <- data.table(Aggregation = setdiff(
    clusterSummaryTable[, unique(Aggregation)], "Raw"))
  rawResult <- clusterSummaryTable[Aggregation == "Raw"][, -"Aggregation"]
  aggRawTable <- aggRawTable[, as.list(rawResult), by = Aggregation]
  # Put Aggregation column back to its original position
  setcolorder(aggRawTable, c(2:4, 1, 5:ncol(clusterSummaryTable)))
  result <- rbind(aggRawTable, clusterSummaryTable[Aggregation != "Raw"])
  return(formatAndSortClusterSummaryTable(result))
}

# Inverts "addRawResultsToClusterSummaryTable()"
revertRawResultsInClusterSummaryTable <- function(clusterSummaryTable) {
  rawLevel <- clusterSummaryTable[, max(Level)]
  aggRawTable <- clusterSummaryTable[Level == rawLevel]
  # We can take any of the aggregation, because all are same (copied)
  aggRawTable <- aggRawTable[Aggregation == Aggregation[1]]
  aggRawTable[, Aggregation := "Raw"]
  result <- rbind(aggRawTable, clusterSummaryTable[Level != rawLevel])
  return(formatAndSortClusterSummaryTable(result))
}

equalOrBothNA <- function(x,y) {
  return(x == y | (is.na(x)) & is.na(y))
}

# Checks if two data tables (in our context: cluster summary tables, but the
# routine is generic) have the same content.
checkTableIdentity <- function(clusterSummaryTable1, clusterSummaryTable2) {
  if (nrow(clusterSummaryTable1) != nrow(clusterSummaryTable2)) {
    cat("Different number of rows.\n")
    return(FALSE)
  }
  if (ncol(clusterSummaryTable1) != ncol(clusterSummaryTable2)) {
    cat("Different number of columns.\n")
    return(FALSE)
  }
  if (any(colnames(clusterSummaryTable1) != colnames(clusterSummaryTable2))) {
    cat("Different column names.\n")
    return(FALSE)
  }
  columnNames <- colnames(clusterSummaryTable1)
  differentValues <- sapply(columnNames, function(colName) {
    return(!all(equalOrBothNA(clusterSummaryTable1[, get(colName)],
                 clusterSummaryTable2[, get(colName)])))
  })
  errorCount <- sum(differentValues)
  if (errorCount > 0) {
    if (errorCount > 10) {
      cat("The following columns have different values: ",
          paste0(columnNames[differentValues], collapse = ", "), "\n.", sep = "")
    } else {
      cat(errorCount, " out of ", ncol(clusterSummaryTable1), " columns have ",
          "different values.\n", sep = "")
    }
    return(FALSE)
  }
  return(TRUE)
}

# Converts the columns of a clustering results (assignments, evaluation) table
# and sorts them
formatAndSortClusterSummaryTable <- function(clusterSummaryTable) {
  result <- copy(clusterSummaryTable)
  # Column conversion - as.character first so levels are sorted alphabetically
  result[, Dataset := as.factor(as.character(Dataset))]
  result[, Dissimilarity := as.factor(as.character(Dissimilarity))]
  result[, Clustering := as.factor(as.character(Clustering))]
  result[, Aggregation := as.factor(as.character(Aggregation))]
  result[, Level := as.numeric(Level)]
  setorder(result, Dataset, Dissimilarity, Clustering, Aggregation, -Level)
  return(result)
}

# Converts the columns of a clustering results (assignments, evaluation) table,
# using character instead of factor for descriptive columns so joins with other
# tables (with different levels) are possible
formatClusterSummaryTableForJoin <- function(clusterSummaryTable) {
  result <- copy(clusterSummaryTable)
  result[, Dataset := as.character(Dataset)]
  result[, Dissimilarity := as.character(Dissimilarity)]
  result[, Clustering := as.character(Clustering)]
  result[, Aggregation := as.character(Aggregation)]
  result[, Level := as.numeric(Level)]
  setorder(result, Dataset, Dissimilarity, Clustering, Aggregation, -Level)
  return(result)
}

# Takes a list of paths to dissimilarity matrices, containing one file with the
# aggregation type "Raw" for each Clustering, Dissimilarity, Dataset and replaces
# them with multiple matrices (files), one for each aggregation type (so we have
# base values for each aggregation and other evaluation methods can work with it)
addRawResultsToDissMatrices <- function(dissMatrixPathList) {
  # naming scheme depends on analyzeAggregationImpactSingleDis()
  for (dissMatrixPath in dissMatrixPathList) {
    dissMatrixFiles <- list.files(path = dissMatrixPath,
        pattern = "^DissMatrix_.*\\.rds$", full.names = TRUE)
    dissMatrixNames <- gsub("\\.rds", "", basename(dissMatrixFiles))
    aggregationNames <- unique(sapply(strsplit(dissMatrixNames, "_"),
                                      function(x) x[[4]]))
    aggregationNames <- setdiff(aggregationNames, "Raw")
    rawDissMatrixFiles <- grep("Raw", dissMatrixFiles, value = TRUE)
    for (rawDissMatrixFile in rawDissMatrixFiles) {
      for (aggName in aggregationNames) {
        file.copy(from = rawDissMatrixFile, to = gsub(pattern = "Raw",
            replacement = aggName, x = rawDissMatrixFile), overwrite = TRUE)
      }
      file.remove(rawDissMatrixFile)
    }
  }
  return(TRUE)
}

### General evaluation ###

# Checks if the design choice columns/experimental settings are the same in two
# evaluation tables. Pre-requisite for calculating correlations.
checkDesignColumnIdentity <- function(clusterSummaryTable1, clusterSummaryTable2) {
  if (clusterSummaryTable1[, .N] != clusterSummaryTable2[, .N]) {
    cat("Different number of rows")
    return(FALSE)
  }
  if (any(clusterSummaryTable1$Dataset != clusterSummaryTable2$Dataset)) {
    cat("Different datasets or dataset order.")
    return(FALSE)
  }
  if (any(clusterSummaryTable1$Dissimilarity != clusterSummaryTable2$Dissimilarity)) {
    cat("Different dissimilarities or dissimilarity order.")
    return(FALSE)
  }
  if (any(clusterSummaryTable1$Clustering != clusterSummaryTable2$Clustering)) {
    cat("Different clusterings or clustering order.")
    return(FALSE)
  }
  if (any(clusterSummaryTable1$Aggregation != clusterSummaryTable2$Aggregation)) {
    cat("Different aggregations or aggregation order.")
    return(FALSE)
  }
  if (any(clusterSummaryTable1$Level != clusterSummaryTable2$Level)) {
    cat("Different levels or level order.")
    return(FALSE)
  }
  return(TRUE)
}

# For each experimental setting (except those in the column "categoryName"), find
# the best setting of the column "categoryName"; attention: if there is only one
# value for a certain setting, is is not included in the result, as we cannot
# make a comparison
createBestExperimentSettingsTable <- function(clusterSummaryTable, categoryName,
    maximize = TRUE, na.rm = FALSE) {
  clusterSummaryTable <- copy(clusterSummaryTable)
  groupCols <- setdiff(colnames(clusterSummaryTable), c(categoryName, "Value"))
  if (na.rm) {
    clusterSummaryTable <- clusterSummaryTable[!is.na(Value)]
  }
  if (categoryName == "Level") { # Base Level should be considered for all aggregations
    clusterSummaryTable <- addRawResultsToClusterSummaryTable(clusterSummaryTable)
  }
  # The experimental settings might be no true cross-product, but for some
  # settings there might be only one value for the category we are interested in
  # (e.g. Aggregation="Raw" would always win on the highest Level); we don't want
  # to count wins there, as winning without competition is kind of cheap (and misleading)
  multiValueSettings <- clusterSummaryTable[, .N > 1, by = groupCols][V1 == TRUE][, -"V1"]
  clusterSummaryTable <- merge(clusterSummaryTable, multiValueSettings,
                               by = colnames(multiValueSettings))
  if (maximize) {
    # not using which.max because there might be several maxima
    winnerTable <- clusterSummaryTable[, get(categoryName)[Value == max(Value)],
                                       by = groupCols]
  } else {
    winnerTable <- clusterSummaryTable[, get(categoryName)[Value == min(Value)],
                                       by = groupCols]
  }
  setnames(winnerTable, "V1", categoryName)
  return(winnerTable)
}

# Counts which setting (dissimilarity, aggregation etc.) yields the highest or
# lowest value in an evaluation table (grouping by all remaining columns, i.e.,
# experimental settings + index + reference)
countBestExperimentSetting <- function(clusterSummaryTable, categoryName, maximize = TRUE,
    relativeCount = FALSE, na.rm = FALSE) {
  groupCols <- setdiff(colnames(clusterSummaryTable), c(categoryName, "Value"))
  winnerTable <- createBestExperimentSettingsTable(
    clusterSummaryTable = clusterSummaryTable, categoryName = categoryName,
    maximize = maximize, na.rm = na.rm)
  if (relativeCount) {
    settingsCount <- nrow(winnerTable[, 1, by = groupCols])
    result <- winnerTable[, .(N = .N / settingsCount), by = categoryName]
  } else {
    result <- winnerTable[, .N, by = categoryName]
  }
  setorder(result, -N)
  return(result)
}

# Calculates the Mean Level Deviation (MLD) expressing the volatility over the
# aggregation levels; averages by the columns provided in "categories"; if the
# index you are interested in is not in the column "Value", provide the param
# "valueCol"
createMLDTable <- function(evaluationTable, categories, valueCol = "Value") {
  if (evaluationTable[Aggregation == "Raw", .N] > 1) {
    evaluationTable <- addRawResultsToClusterSummaryTable(evaluationTable)
  }
  if ("Index" %in% colnames(evaluationTable) && !("Index" %in% categories) &&
      evaluationTable[, uniqueN(Index)] > 1) {
    warning("You are also averaging over indices. Are you sure?")
  }
  if ("Reference" %in% colnames(evaluationTable) && !("Reference" %in% categories) &&
      evaluationTable[, uniqueN(Reference)] > 1) {
    warning("You are also averaging over references. Are you sure?")
  }
  groupCols <- intersect(colnames(evaluationTable), c(categories, "Dataset",
      "Aggregation", "Dissimilarity", "Clustering", "Reference", "Index"))
  return(evaluationTable[, .(Stddev = sd(get(valueCol))), by = groupCols][,
      .(MLD = mean(Stddev)), by = categories][order(MLD)])
}

### External validation (relative) ###

# Calculates several external validity indices for two cluster assignments,
# returning a data.table containing the index names and values
calcExternalCVIs <- function(clusterAssignments1, clusterAssignments2) {
  pairCVIParams <- pairCVIParameters_fast(clusterAssignments1, clusterAssignments2)
  return(data.table(
    Index = c("Rand", "Fowlkes-Mallows", "Phi", "inv.vanDongen", "NMI"),
    Value = c(
      randIndex_fast(pairCVIParams, normalize = TRUE),
      fowlkesMallows_fast(pairCVIParams, normalize = TRUE),
      phi_fast(pairCVIParams),
      vanDongen_fast(clusterAssignments1, clusterAssignments2, normalizeAndInvert = TRUE),
      VI_fast(clusterAssignments1, clusterAssignments2, normalizeAndInvert = TRUE)
    )
  ))
}

# Calculates several external validity indices for a data.table containing
# multiple aggregation levels (columns "Level" and assignments, first row
# considered as base aggregation/raw data); returns a data.table with the
# columns "Level", "Index", "Value" and "Reference" ("B", "P" or "M")
calcExternalCVIsForAggregation <- function(aggClusterAssignmentsTable) {
  setorder(aggClusterAssignmentsTable, -Level)
  # Clustering compared to base aggregation level (might be raw data)
  baseAggAssignment <- as.numeric(aggClusterAssignmentsTable[1, -"Level"])
  baseAggCVIs <- aggClusterAssignmentsTable[, calcExternalCVIs(baseAggAssignment,
      as.numeric(.SD)), by = Level]
  baseAggCVIs[, Reference := "B"]
  # Clustering compared to previous aggregation level
  if (nrow(aggClusterAssignmentsTable) > 1) {
    assMatrixPrev <- as.matrix(aggClusterAssignmentsTable[-.N, -"Level"])
    assMatrixCur <- as.matrix(aggClusterAssignmentsTable[-1, -"Level"])
    prevAggCVIs <- rbindlist(lapply(1:aggClusterAssignmentsTable[, .N - 1], function(i) {
      curTable <- calcExternalCVIs(assMatrixPrev[i, ], assMatrixCur[i, ])
      curTable[, Level := aggClusterAssignmentsTable[i + 1, Level]]
    }))
    prevAggCVIs <- rbind(
      data.table(Level = aggClusterAssignmentsTable[1, Level],
                 Index = prevAggCVIs[, unique(Index)], Value = NA),
      prevAggCVIs) # base level has no predecessor, therefore we add NAs
  } else {
    warning("Only one level for an aggregation - writing NAs for previous level comparison.")
    prevAggCVIs <- data.table(Level = aggClusterAssignmentsTable[1, Level],
                              Index = baseAggCVIs[, unique(Index)], Value = NA)
  }
  prevAggCVIs[, Reference := "P"]
  # Clustering compared to median aggregation level
  medianAggIndex <- floor((aggClusterAssignmentsTable[, .N] + 1) / 2)
  medianAggAssignment <- as.numeric(aggClusterAssignmentsTable[medianAggIndex, -"Level"])
  medianAggCVIs <- aggClusterAssignmentsTable[, calcExternalCVIs(medianAggAssignment,
      as.numeric(.SD)), by = Level]
  medianAggCVIs[, Reference := "M"]
  return(rbind(baseAggCVIs, prevAggCVIs, medianAggCVIs))
}

# Calculates several external validity indices:
# - takes a data.table containing columns "Dataset", "Clustering", Dissimilarity",
# "Aggregation", "Level" and the clustering assignments (all remaining columns)
# - creates as data.table with the same initial columns and additional columns containing
# the selected CVIs, comparing each level (per dissimilarity and aggregation method) to
# the previous one (P), the initial one (B) and the median level (M)
createExternalCVITable <- function(clusterAssignmentsTable) {
  nonAssignmentAndLevelColumns <- setdiff(colnames(clusterAssignmentsTable),
      c("Level", grep(pattern = "^A[0-9]+", x = colnames(clusterAssignmentsTable),
                      value = TRUE)))
  result <- clusterAssignmentsTable[, calcExternalCVIsForAggregation(.SD),
      by = nonAssignmentAndLevelColumns]
  # Swap "Reference" with "Value" column ("Value" should come last)
  setcolorder(result, c(seq_len(ncol(result) - 2), ncol(result), ncol(result) - 1))
  result[, Index := as.factor(Index)]
  result[, Reference := as.factor(Reference)]
  return(formatAndSortClusterSummaryTable(result))
}

# Prints the number of unexpected NAs in an external CVI table.
# (simplified check, might not detect all errors, but easy to implement)
checkExternalCVINAs <- function(externalCVITable, clusterCountTable) {
  naTable <- externalCVITable[is.na(Value)]
  # Highest level has no previous level
  highestLevel <- externalCVITable[, max(Level)]
  naTable <- naTable[!(Level == highestLevel & Reference == "P")]
  # All normalized indices are NaN if both clusterings consist of only one cluster
  # (we make a more simple check here which might prune more)
  naTable <- merge(naTable, clusterCountTable[Value != 1, -"Value"])
  # Phi index also NaN if one of the clusterings consists of only one cluster
  cat(naTable[Index != "Phi", .N], "unexpected NAs for external validity indices.\n")
  return(TRUE)
}

### External validation (ground truth) ###

# Calculates several external validity indices (quite similar to calcExternalCVIs())
# for a cluster assignment and a ground truth vector, returning a data.table
# containing the index names and values
calcGroundTruthIndices <- function(clusterAssignments, groundTruth) {
  pairCVIParams <- pairCVIParameters_fast(clusterAssignments, groundTruth)
  return(data.table(
    Index = c("Uniformity.clusters", "Uniformity.classes", "Purity.clusters",
              "Purity.classes", "Rand", "Fowlkes-Mallows", "Phi", "inv.vanDongen", "NMI"),
    Value = c(
      conditionalEntropy_fast(clusterAssignments, groundTruth, normalizeAndInvert = TRUE),
      conditionalEntropy_fast(groundTruth, clusterAssignments, normalizeAndInvert = TRUE),
      purity_fast(clusterAssignments, groundTruth),
      purity_fast(groundTruth, clusterAssignments),
      randIndex_fast(pairCVIParams, normalize = TRUE),
      fowlkesMallows_fast(pairCVIParams, normalize = TRUE),
      phi_fast(pairCVIParams),
      vanDongen_fast(clusterAssignments, groundTruth, normalizeAndInvert = TRUE),
      VI_fast(clusterAssignments, groundTruth, normalizeAndInvert = TRUE)
    )
  ))
}

# Prints the number of unexpected NAs in an external CVI table.
# (simplified check, might not detect all errors, but easy to implement)
checkGroundTruthNAs <- function(groundTruthTable) {
  # similar to external CVIs, but (usually) no danger that the ground truth
  # consists of only one class
  cat(groundTruthTable[Index != "Phi", sum(is.na(Value))],
      "unexpected NAs for ground truth indices.\n")
  return(TRUE)
}

# Calculates external validity indices comparing cluster assignments to the ground truth:
# - takes a data.table containing columns "Dataset", "Clustering", Dissimilarity",
# "Aggregation", "Level" and the clustering assignments (all remaining columns)
# - takes an integer vector with values in 1:k', representing a ground truth
# classification of the assigned time series (e.g. day of week, machine etc.)
# - creates as data.table with the same initial columns and additional columns containing
# the external CVIs, comparing each assignment to the ground truth
createGroundTruthIndexTable <- function(clusterAssignmentsTable, groundTruth) {
  groundTruthCategories <- sort(unique(groundTruth))
  if (length(groundTruthCategories) != max(groundTruthCategories) ||
      any(groundTruthCategories != 1:max(groundTruthCategories))) {
    stop("The groundTruth has to be an integer vector with values from 1 to
         max(groundTruth) without gaps.")
  }
  return(formatAndSortClusterSummaryTable(
    clusterAssignmentsTable[, calcGroundTruthIndices(as.numeric(.SD), groundTruth),
        by = setdiff(colnames(clusterAssignmentsTable),
                     grep("^A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE))])
  )
}

### Interval validation ###

# Calculates several internal validity indices for a vector of (integer) cluster
# assignments and a (symmetric) dissimilarity/distance matrix; returns a table
calcInternalCVIs <- function(clusterAssignments, dissMatrix) {
  clusterAssignments <- as.integer(clusterAssignments)
  if (all(clusterAssignments == 1)) { # one cluster contains all objects
    values <- c(NA, 1, NA, NA)
  } else {
    clDissList <- cls.scatt.diss.mx(diss.mx = dissMatrix, clust = clusterAssignments)
    values <- c(
      avgSilhouette(clusterAssignments, dissMatrix),
      normConnectivity(clusterAssignments, dissMatrix, k = 10),
      generalizedDunn_fast(clDissList$intercls.average, clDissList$intracls.average),
      generalizedDB_fast(clDissList$intercls.average, clDissList$intracls.average),
      iGeneralizedDB_fast(clDissList$intercls.average, clDissList$intracls.average)
    )
  }
  return(data.table(
    Index = c("Silhouette", "inv.Connectivity", "Gen.Dunn", "Gen.Davies.Bouldin",
              "inv.Gen.Davies.Bouldin"),
    Value = values
  ))
}

# Calculates several internal validity indices:
# - takes a data.table containing columns "Dataset", "Clustering", Dissimilarity",
# "Aggregation", "Level" and the clustering assignments (all remaining columns)
# - takes a list of paths to directories where dissimilarity matrices are saved as
# *.rds files
# - creates a data.table with the same initial columns and additional columns
# "Index", "Value", "Reference" containing the selected CVIs (C = current CVI =
# CVI for that assignment and diss matrix, B = base CVI = CVI for that assignment,
# but using diss matrix of first aggregation step [usually raw data])
createInternalCVITable <- function(clusterAssignmentsTable, dissMatrixPathList) {
  nonAssignmentColumns <- setdiff(colnames(clusterAssignmentsTable),
      grep(pattern = "^A[0-9]+", x = colnames(clusterAssignmentsTable), value = TRUE))
  dissMatrixPaths <- unlist(lapply(dissMatrixPathList, list.files,
      pattern = "^DissMatrix_.*\\.rds$", full.names = TRUE))
  dissMatrixNames <- gsub("\\.rds", "", basename(dissMatrixPaths))
  # CVIs using assignments and their corresponding dissimilarity matrices
  computingCluster <- makeCluster(detectCores())
  registerDoParallel(computingCluster)
  currentAggCVIs <- rbindlist(foreach(i = 1:length(dissMatrixPaths),
      .packages = c("cluster", "clv", "data.table", "FastTSDistances"),
      .export = c("avgSilhouette", "calcInternalCVIs", "normConnectivity")) %dopar% {
    # assumed naming convection (see analyzeAggregationImpactSingleDis()):
    # <<irrelevant>>_<<Dataset>>_<<Dissimilarity>>_<<Aggregation>>_<<Level>>
    matrixNameParts <- strsplit(dissMatrixNames[i], "_")[[1]]
    relevantAssignments <- clusterAssignmentsTable[Dataset == matrixNameParts[2] &
        Dissimilarity == matrixNameParts[3] & Aggregation == matrixNameParts[4] &
        Level == matrixNameParts[5]]
    if (nrow(relevantAssignments) == 0) {
      stop(paste0("No suitable entry in clusterAssignmentsTable for dissimilarity ",
                  "matrix ", dissMatrixNames[i]))
    }
    dissMatrix <- readRDS(dissMatrixPaths[i])
    return(relevantAssignments[, calcInternalCVIs(as.numeric(.SD), dissMatrix),
                               by = nonAssignmentColumns])
  })
  stopCluster(computingCluster)
  # Quick sanity check: for each clus result there should be multiple CVIs
  if (nrow(currentAggCVIs) %% nrow(clusterAssignmentsTable) != 0) {
    stop("Some cluster assignments seem to be related to no dissimilarity matrix.")
  }
  currentAggCVIs[, Reference := "C"] # new column, put at 2nd-to-last position
  setcolorder(currentAggCVIs, c(seq_len(ncol(currentAggCVIs) - 2),
                                ncol(currentAggCVIs), ncol(currentAggCVIs) - 1))
  # CVIs using assignments and the base aggregation dissimilarity matrices
  baseAggregationLevel <- clusterAssignmentsTable[, max(Level)]
  baseDissMatrixPaths <- grep(paste0("_", baseAggregationLevel, "\\.rds$"),
                                    dissMatrixPaths, value = TRUE)
  baseDissMatrixNames <- gsub("\\.rds", "", basename(baseDissMatrixPaths))
  computingCluster <- makeCluster(detectCores())
  registerDoParallel(computingCluster)
  baseAggCVIs <- rbindlist(foreach(i = 1:length(baseDissMatrixPaths),
      .packages = c("cluster", "clv", "data.table", "FastTSDistances"),
      .export = c("avgSilhouette", "calcInternalCVIs", "normConnectivity")) %dopar% {
    matrixNameParts <- strsplit(baseDissMatrixNames[i], "_")[[1]]
    # difference to code above: "Level" not more relevant for sub-setting
    relevantAssignments <- clusterAssignmentsTable[Dataset == matrixNameParts[2] &
        Dissimilarity == matrixNameParts[3]]
    # If there is a special "Raw" pseudo-aggregation, use its diss matrix as base,
    # else the highest aggregation level should exist for each aggregation type separately
    if (matrixNameParts[4] != "Raw") {
      relevantAssignments <- relevantAssignments[Aggregation == matrixNameParts[4]]
    }
    dissMatrix <- readRDS(baseDissMatrixPaths[i])
    return(relevantAssignments[, calcInternalCVIs(as.numeric(.SD), dissMatrix),
                               by = nonAssignmentColumns])
  })
  stopCluster(computingCluster)
  baseAggCVIs[, Reference := "B"]
  setcolorder(baseAggCVIs, c(seq_len(ncol(currentAggCVIs) - 2),
                             ncol(currentAggCVIs), ncol(currentAggCVIs) - 1))
  # Combine
  result <- rbind(currentAggCVIs, baseAggCVIs)
  result[, Index := as.factor(Index)]
  result[, Reference := as.factor(Reference)]
  return(formatAndSortClusterSummaryTable(result))
}

# Checks if NA values in an internalCVITable can be explained and prints unusual
# values; always returns TRUE
checkInternalCVINAs <- function(internalCVITable, clusterAssignmentsTable) {
  clusterCountTable <- createClusterCountTable(clusterAssignmentsTable)
  allSingletonClusterTable <- createAllSingletonClustersTable(clusterAssignmentsTable)
  # Silhouette: NA if there is only one cluster or only singleton clusters
  silNATable <- internalCVITable[Index == "Silhouette" & is.na(Value)]
  silNATable <- merge(silNATable, clusterCountTable[Value != 1, -"Value"])
  silNATable <- merge(silNATable, allSingletonClusterTable[Value == FALSE, -"Value"])
  cat(nrow(silNATable), "unexpected NAs for silhouette index.\n")
  # Connectivity: no NAs expected
  cat(internalCVITable[Index == "inv.Connectivity", sum(is.na(Value))], " unexpected ",
      "NAs for inverted connectivity index.\n", sep = "")
  # Generalized Dunn: NA if only one cluster (inter-clus diss not available)
  # (can also be Inf if all intra-cluster diss values are 0, e.g. for EDR and short series)
  dunnNATable <- internalCVITable[Index == "Gen.Dunn" & is.na(Value)]
  dunnNATable <- merge(dunnNATable, clusterCountTable[Value != 1, -"Value"])
  cat(nrow(dunnNATable), "unexpected NA for Generalized Dunn index.\n")
  # Generalized Davies-Bouldin: NA if only one cluster or we have clusters
  # with inter-cluster diss == 0 and intra-clus diss == 0 (affinity propagation
  # and EDR dissimilarity prone to this behaviour, but might also happen for
  # comparisons on base time series combined with clustering result on agg series)
  # (can also be Inf if there is a inter-clus diss == 0 with intra-cluss diss != 0)
  dbNATable <- internalCVITable[Index == "Gen.Davies.Bouldin" & is.na(Value)]
  dbNATable <- merge(dbNATable, clusterCountTable[Value != 1, -"Value"])
  cat(nrow(dbNATable), " unexpected NAs for Generalized Davies-Bouldin ",
      "index, including ", dbNATable[Clustering == "AffProp", .N], " for affinity ",
      "propagation, ", dbNATable[Reference == "B", .N], " for base reference and ",
      dbNATable[Dissimilarity == "EDR", .N], " for EDR.\n", sep = "")
  # Inverted Generalized Davies-Bouldin: NA if only one cluster or we have clusters
  # with inter-cluster diss == 0 and intra-clus diss == 0 (affinity propagation
  # and EDR dissimilarity prone to this behaviour, but might also happen for
  # comparisons on base time series)
  # (can also be Inf if there are at least two clusters with intra-diss == 0)
  dbNATable <- internalCVITable[Index == "inv.Gen.Davies.Bouldin" & is.na(Value)]
  dbNATable <- merge(dbNATable, clusterCountTable[Value != 1, -"Value"])
  cat(nrow(dbNATable), " unexpected NAs for Inverted Generalized Davies-Bouldin ",
      "index, including ", dbNATable[Clustering == "AffProp", .N], " for affinity ",
      "propagation, ", dbNATable[Reference == "B", .N], " for base reference and ",
      dbNATable[Dissimilarity == "EDR", .N], " for EDR.\n", sep = "")
  return(TRUE)
}

# Takes a data.table containing validity indices as absolute numbers and returns a
# table of the same length containing the differences between CVIs of consectutive
# aggregation levels; if "absChange", the absolute value of the differences will
# be used
transformToCVIChangeTable <- function(clusteringEvaluationTable, absChange = FALSE) {
  result <- copy(clusteringEvaluationTable)
  if (absChange) {
    result[, Value := c(NA, abs(diff(Value))),
           by = setdiff(colnames(result), c("Value", "Level"))]
  } else {
    result[, Value := c(NA, diff(Value)),
           by = setdiff(colnames(result), c("Value", "Level"))]
  }
  return(result)
}

### Entropy ###

# Calculates the shannon entropy of a vector x, considering each distinct value
# as a different class. If "normalize", then the value is divided by the maximum
# possible entropy, resulting in a [0,1] value (or NaN if there is only one class)
entropy <- function(x, normalize = FALSE) {
  relFreq <- table(x) / length(x)
  relFreq <- relFreq[relFreq > 0] # factors with unused levels
  entropy <- -sum(log2(relFreq) * relFreq)
  if (normalize) {
    return(entropy / log2(uniqueN(x)))
  } else {
    return(entropy)
  }
}

# Calculates the entropy of each clustering, returning a data.table with the same
# number of rows, but one "Value" column with the entropy instead of the assignment
# columns
createClusteringEntropyTable <- function(clusterAssignmentsTable, normalize = FALSE) {
  return(formatAndSortClusterSummaryTable(
    clusterAssignmentsTable[, .(Value = clusterEntropy_fast(as.numeric(.SD), normalize)),
        by = setdiff(colnames(clusterAssignmentsTable),
                     grep("^A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE))])
  )
}

# Prints the number of unexpected NAs in an entropy table.
checkEntropyNAs <- function(entropyTable) {
  # no NAs expected (regardless if normalized or not)
  cat(entropyTable[, sum(is.na(Value))], "unexpected NAs for entropy.\n")
  return(TRUE)
}

# Table which has the usual meta-columns + a column with the number of clusters
createClusterCountTable <- function(clusterAssignmentsTable) {
  return(formatAndSortClusterSummaryTable(
    clusterAssignmentsTable[, .(Value = uniqueN(as.numeric(.SD))),
        by = setdiff(colnames(clusterAssignmentsTable),
                     grep("^A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE))])
  )
}

# Table which has the usual meta-columns + a column with the size of the biggest cluster
createMaxClusterSizeTable <- function(clusterAssignmentsTable) {
  return(formatAndSortClusterSummaryTable(
    clusterAssignmentsTable[, .(Value = max(tabulate(as.numeric(.SD)))),
        by = setdiff(colnames(clusterAssignmentsTable),
                     grep("^A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE))])
  )
}

# Table which has the usual meta-columns + a logical column which states if the
# clustering contains any singleton clusters or not
createAnySingletonClustersTable <- function(clusterAssignmentsTable) {
  return(clusterAssignmentsTable[, .(Value = any(tabulate(as.numeric(.SD)) == 1)),
      by = setdiff(colnames(clusterAssignmentsTable),
                   grep("^A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE))])
}

# Table which has the usual meta-columns + a logical column which states if the
# clustering contains only singleton clusters
createAllSingletonClustersTable <- function(clusterAssignmentsTable) {
  return(clusterAssignmentsTable[, .(Value = all(tabulate(as.numeric(.SD)) == 1)),
      by = setdiff(colnames(clusterAssignmentsTable),
                   grep("^A[0-9]+", colnames(clusterAssignmentsTable), value = TRUE))])
}

### Forecasting ###

# Calculates the mean absolute error between a forecasted time series and the
# actual time series (can be univariate or multivariate).
calcMAE <- function(forecast, actual) {
  return(sum(abs(actual - forecast)) / length(forecast))
}


# Calculates the mean absolute percentage error between a forecasted time series
# and the actual time series (can be univariate or multivariate).
calcMAPE <- function(forecast, actual) {
  return(100 / length(forecast) * sum(abs((actual - forecast) / actual)))
}

# Takes an already clustered time series list "oldDatasetPrev" with the cluster
# assignments vector (same length) "oldDataAssignments" and the corresponding
# succeeding time series "oldDatasetActual",
# a test dataset of series to forecast "newDatasetActual" and their predecessors
# (same length) "newDatasetPrev" as list of time series with the assignment vector
# (same length) "newDataAssignments"
# and a dissimilarity matrix of size length(newDataset) * length(oldDataset).
#
# Makes instance-based forecasts for all new time series and calculates the
# forecast error in four categories/forecasting approaches:
# 1) based on all series in "oldDatasetActual", unweighted mean
# 2) based on all series in "oldDatasetActual", weighted with similarities of their
# predecessors "oldDatasetPrev" to the "newDatasetPrev"
# 3) based on series in "oldDatasetActual" belonging to same cluster as the new
# series, unweighted mean
# 4) based on series in "oldDatasetActual" belonging to same cluster as the new
# series, weighted with similarities
# 5) based on the series in "oldDatasetActual" whose predecessor in "oldDatasetPrev"
# is closest to new series in "newDatasetPrev" (1NN forecast)
# 6) 10NN forecast, unweighted mean
# 7) 10NN forecast, weighted with similarities
# 8) based on the series in "newDatasetPrev" which comes before the forecast
# interval (simply using it for the next internal, naive forecast)
#
# The parameter "dissWeightConversion" determines how dissimilarities are turned
# into weights.
calcForecastError <- function(newDatasetPrev, newDatasetActual, oldDatasetPrev,
    oldDatasetActual, newDataAssignments, oldDataAssignments, dissMatrix,
    dissWeightConversion = "minmax") {
  # Error based on all training series, unweighted mean
  allMeanSeries <- averageTimeSeriesList(tsList = oldDatasetActual,
      dissimilarities = rep(1, length(oldDatasetActual)), dissWeightConversion)
  allMeanError <- mean(sapply(1:length(newDatasetActual), function(i) {
    return(calcMAE(forecast = allMeanSeries, actual = newDatasetActual[[i]]))
  }))
  # Error based on all training time series, weighted with dissimilarity
  allWeightedError <- mean(sapply(1:length(newDatasetActual), function(i) {
    # to predict newActual[i], find similarity of oldPrev[..] with newPrev[i] and use
    # it to weight oldActual[..] as forecasts (instance-based prediction)
    forecast <- averageTimeSeriesList(tsList = oldDatasetActual,
        dissimilarities = dissMatrix[i, ], dissWeightConversion)
    return(calcMAE(forecast = forecast, actual = newDatasetActual[[i]]))
  }))
  # Error based on training time series in same cluster, unweighted mean
  clusMeanError <- mean(sapply(1:length(newDatasetActual), function(i) {
    sameClusterIdx <- which(oldDataAssignments == newDataAssignments[i])
    forecast <- averageTimeSeriesList(tsList = oldDatasetActual[sameClusterIdx],
        dissimilarities = rep(1, length(sameClusterIdx)), dissWeightConversion)
    return(calcMAE(forecast = forecast, actual = newDatasetActual[[i]]))
  }))
  # Error based on training time series in same cluster, weighted with dissimilarity
  clusWeightedError <- mean(sapply(1:length(newDatasetActual), function(i) {
    sameClusterIdx <- which(oldDataAssignments == newDataAssignments[i])
    forecast <- averageTimeSeriesList(tsList = oldDatasetActual[sameClusterIdx],
        dissimilarities = dissMatrix[i, sameClusterIdx], dissWeightConversion)
    return(calcMAE(forecast = forecast, actual = newDatasetActual[[i]]))
  }))
  # Error based on closest training time series
  nn1Error <- mean(sapply(1:length(newDatasetActual), function(i) {
    nnIdx <- which.min(dissMatrix[i, ])
    forecast <- oldDatasetActual[[nnIdx]]
    return(calcMAE(forecast = forecast, actual = newDatasetActual[[i]]))
  }))
  # Error based on 10 closest training time series, unweighted mean
  nn10MeanError <- mean(sapply(1:length(newDatasetActual), function(i) {
    nnIdx <- order(dissMatrix[i, ])[1:10]
    forecast <- averageTimeSeriesList(tsList = oldDatasetActual[nnIdx],
        dissimilarities = dissMatrix[i, nnIdx], dissWeightConversion)
    return(calcMAE(forecast = forecast, actual = newDatasetActual[[i]]))
  }))
  # Error based on 10 closest training time series, weighted with dissimilarity
  nn10WeightedError <- mean(sapply(1:length(newDatasetActual), function(i) {
    nnIdx <- order(dissMatrix[i, ])[1:10]
    forecast <- averageTimeSeriesList(tsList = oldDatasetActual[nnIdx],
        dissimilarities = rep(1, 10), dissWeightConversion)
    return(calcMAE(forecast = forecast, actual = newDatasetActual[[i]]))
  }))
  # Error using the time series itself as prediction for the next interval
  naiveError <- mean(sapply(1:length(newDatasetActual), function(i) {
    return(calcMAE(forecast = newDatasetPrev[[i]], actual = newDatasetActual[[i]]))
  }))
  return(data.table(
    Index = rep("MAE", 8),
    Reference = c("all.mean", "all.weighted", "clus.mean", "clus.weighted",
                  "1NN", "10NN.mean", "10NN.weighted", "naive"),
    Value = c(allMeanError, allWeightedError, clusMeanError, clusWeightedError,
              nn1Error, nn10MeanError, nn10WeightedError, naiveError)
  ))
}

# Graphical routine to create forecasts for all members of a test set based on
# a training set, cluster assignments and dissimilarities (parameters and
# computation similar to "calcForecastError()"). Can handle univariate as well
# as multivariate time series (but all time series should have the same number
# of attributes and same attribute names).
# "withBackground" determines if all time series which are combined in the
# clustering approach should be plotted in gray as background
plotForecasts <- function(newDatasetPrev, newDatasetActual, oldDatasetPrev,
    oldDatasetActual, newDataAssignments, oldDataAssignments, dissMatrix,
    plotDir = "./", dissWeightConversion = "minmax", attributeName = NULL,
    withBackground = TRUE, progressBar = TRUE) {
  if (!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = TRUE)
  }
  if (is.matrix(newDatasetActual[[1]])) {# multivariate time series -> plots for each attribute
    for (j in 1:ncol(newDatasetActual[[1]])) {
      attributeName <- colnames(newDatasetActual[[1]])[j]
      if (is.null(attributeName)) {
        attributeName <- paste0("Attribute_", j)
      }
      plotForecasts(newDatasetPrev = lapply(newDatasetPrev, function(x) x[, j]),
                    newDatasetActual = lapply(newDatasetActual, function(x) x[, j]),
                    oldDatasetPrev = lapply(oldDatasetPrev, function(x) x[, j]),
                    oldDatasetActual = lapply(oldDatasetActual, function(x) x[, j]),
                    newDataAssignments, oldDataAssignments,
                    dissMatrix, plotDir, dissWeightConversion, attributeName,
                    withBackground)
    }
  } else {
    itemNrs <- seq_len(length(newDatasetActual[[1]]))
    forecastCount <- length(newDatasetActual)
    computingCluster <- makeCluster(detectCores())
    registerDoParallel(computingCluster)
    foreach(i = 1:forecastCount, .packages = c("data.table", "FastTSDistances",
                                               "ggplot2", "tcltk")) %dopar% {
      if (progressBar) {
        if (!exists("progBar")) {
          progBar <- tkProgressBar(max = forecastCount, title = "Plot progress")
        }
        setTkProgressBar(progBar, value = i, label = paste0(ifelse(is.null(attributeName),
            "", paste0("Attribute: ", attributeName, " - ")), "Creating plot ", i,
             "/", forecastCount))
      }
      actualData <- data.table(Item = itemNrs, Value = newDatasetActual[[i]],
                               Type = "Actual")
      allMeanData <- data.table(Item = itemNrs, Value = averageTimeSeriesList(
        tsList = oldDatasetActual, dissimilarities = rep(1, length(itemNrs))),
        Type = "All.mean")
      allWeightedData <- data.table(Item = itemNrs, Value = averageTimeSeriesList(
        tsList = oldDatasetActual, dissimilarities = dissMatrix[i, ],
        dissWeightConversion = dissWeightConversion), Type = "All.weighted")
      sameClusterIdx <- which(oldDataAssignments == newDataAssignments[i])
      clusterMeanData <- data.table(Item = itemNrs, Value = averageTimeSeriesList(
        tsList = oldDatasetActual[sameClusterIdx],
        dissimilarities = rep(1, length(sameClusterIdx))), Type = "Cluster.mean")
      clusterWeightedData <- data.table(Item = itemNrs, Value = averageTimeSeriesList(
        tsList = oldDatasetActual[sameClusterIdx],
        dissimilarities = dissMatrix[i, sameClusterIdx],
        dissWeightConversion = dissWeightConversion), Type = "Cluster.weighted")
      nn1Data <- data.table(Item = itemNrs, Value = oldDatasetActual[[which.min(
        dissMatrix[i, ])]], Type = "1NN")
      nnIdx <- order(dissMatrix[i, ])[1:10]
      nn10MeanData <- data.table(Item = itemNrs, Value = averageTimeSeriesList(
        tsList = oldDatasetActual[nnIdx], dissimilarities = rep(1, length(nnIdx))),
        Type = "10NN.mean")
      nn10WeightedData <- data.table(Item = itemNrs, Value = averageTimeSeriesList(
        tsList = oldDatasetActual[nnIdx], dissimilarities = dissMatrix[i, nnIdx]),
        Type = "10NN.weighted")
      naiveForecastData <- data.table(Item = itemNrs, Value = newDatasetPrev[[i]],
                                      Type = "Naive")
      plotData <- rbind(actualData, allMeanData, allWeightedData, clusterMeanData,
          clusterWeightedData, nn1Data, nn10MeanData, nn10WeightedData, naiveForecastData)
      forecastPlot <- ggplot()
      if (withBackground) {
        forecastPlot <- forecastPlot + lapply(sameClusterIdx, function(idx) {
          geom_line(data = data.table(Item = itemNrs, Value = oldDatasetActual[[idx]]),
                    mapping = aes(x = Item, y = Value), color = "darkgrey")})
      }
      forecastPlot <- forecastPlot +
          geom_line(data = plotData, mapping = aes(x = Item, y = Value, color = Type,
                                                   linetype = Type), size = 1) +
          scale_color_manual("Type", values = c("Actual" = "black", "All.mean" = "blue",
              "All.weighted" = "blue", "Cluster.mean" = "red", "Cluster.weighted" = "red",
              "1NN" = "darkgreen", "10NN.mean" = "brown", "10NN.weighted" = "brown",
              "Naive" = "darkgreen")) +
          scale_linetype_manual("Type", values = c("Actual" = "solid", "All.mean" = "dashed",
              "All.weighted" = "solid", "Cluster.mean" = "dashed", "Cluster.weighted" = "solid",
              "1NN" = "solid", "10NN.mean" = "dashed", "10NN.weighted" = "solid",
              "Naive" = "dashed")) +
          ggtitle(paste0("Forecasts for test series ", i, ifelse(is.null(attributeName),
              "", paste0(" - Attribute: ", attributeName))))
      png(filename = paste0(plotDir, "Forecast_", i, ifelse(is.null(attributeName), "",
          paste0("_Attribute_", attributeName)), ".png"), width = 960, height = 540)
      print(forecastPlot)
      dev.off()
      return(TRUE)
    }
    stopCluster(computingCluster)
  }
  return(TRUE)
}

# Summarizes forecasting errors on two (raw) datasets with and without clustering:
# - takes two lists of list of time series (train/test parts of multiple datasets)
# called "oldDatasetList"/"newDatasetList"
# - takes two lists of lists of two index vectors ("prev" and "actual") which
# determine which time series act as predecessors in the forecasting process and
# which time series act as exemplars for the target to forecast; the outer list
# should be named like the datasets; will be used for sub-setting the time series
# lists (input part, target part in forecasting)
# - takes two data.tables containing columns "Dataset", "Clustering", Dissimilarity",
# "Aggregation", "Level" and the clustering assignments (all remaining columns);
# "oldDataAssignmentsTable" should be the result of "analyzeAggregationImpact()",
# "newDataAssignmentsTable" should be the result of "assignDatasetsToClusters()"
# - takes a vector of directory paths to dissimilarity matrices, describing
# the dissimilarity of the test time series to the training time series (created
# by the assignment routine "assignDatasetsToClusters()")
# - creates a data.table with the same initial columns as the assignment tables
# and additional columns "Index", "Reference", "Value", containing the forecasting
# errors on the datasets (without aggregation) using/not using the cluster
# assignments and dissimilarities for weighting
createBaseForecastErrorTable <- function(newDatasetList, oldDatasetList,
    newDataPrevActList, oldDataPrevActList, newDataAssignmentsTable,
    oldDataAssignmentsTable, dissMatrixPaths, dissWeightConversion = "minmax") {
  # Change column types (should be same in both tables, avoid factors)
  newDataAssignmentsTable <- formatClusterSummaryTableForJoin(newDataAssignmentsTable)
  oldDataAssignmentsTable <- formatClusterSummaryTableForJoin(oldDataAssignmentsTable)
  # Create a joint assignmensTable because old-new pairs are used below
  oldDataAssCols <- grep(pattern = "^A[0-9]+", x = colnames(oldDataAssignmentsTable),
                         value = TRUE)
  newDataCols <- colnames(newDataAssignmentsTable)
  newDataOldAssCols <- grep(pattern = "^A[0-9]+", x = newDataCols, value = TRUE)
  newDataNewAssCols <- gsub("^A", "A\\.", newDataOldAssCols)
  nonAssignmentColumns <- setdiff(newDataCols, newDataOldAssCols)
  # Rename so we can distinguish old and new assignments
  setnames(newDataAssignmentsTable, old = newDataOldAssCols, new = newDataNewAssCols)
  newDataRowCount <- nrow(newDataAssignmentsTable)
  newDataAssignmentsTable <- merge(newDataAssignmentsTable, oldDataAssignmentsTable,
                                   all.x = TRUE)
  if (nrow(newDataAssignmentsTable) != newDataRowCount) {
    stop("There is a Dataset-Clustering-Dissimilarity-Aggregation-Level combination
         in newDataAssignmentsTable which does not exist in oldDataAssignmentsTable.")
  }
  # Find dissimilarity matrix files
  dissMatrixPaths <- list.files(dissMatrixPaths, pattern = "^DissMatrix_.*\\.rds$",
                                full.names = TRUE)
  dissMatrixNames <- gsub("\\.rds", "", basename(dissMatrixPaths))
  dissMatrixNameParts <- strsplit(dissMatrixNames, "_")
  # Forecasts using assignments and their corresponding dissimilarity matrices
  computingCluster <- makeCluster(detectCores())
  registerDoParallel(computingCluster)
  result <- rbindlist(foreach(i = 1:length(dissMatrixPaths),
      .packages = c("data.table", "FastTSDistances"),
      .export = c("calcMAE", "calcForecastError")) %dopar% {
    # assumed naming convection (see analyzeAggregationImpactSingleDis()):
    # <<irrelevant>>_<<Dataset>>_<<Dissimilarity>>_<<Aggregation>>_<<Level>>
    matrixNameParts <- dissMatrixNameParts[[i]]
    datasetName <- matrixNameParts[2]
    relevantAssignments <- newDataAssignmentsTable[Dataset == datasetName &
        Dissimilarity == matrixNameParts[3] & Aggregation == matrixNameParts[4] &
          Level == matrixNameParts[5]]
    if (nrow(relevantAssignments) == 0) {
      stop(paste0("No suitable entry in newDataAssignmentsTable for dissimilarity ",
                  "matrix ", dissMatrixNames[i]))
    }
    dissMatrix <- readRDS(dissMatrixPaths[i])
    # Forecasts are made based on time series and their predecessors
    newPrevIdx <- newDataPrevActList[[datasetName]]$prev
    newActIdx <- newDataPrevActList[[datasetName]]$actual
    oldPrevIdx <- oldDataPrevActList[[datasetName]]$prev
    oldActIdx <- oldDataPrevActList[[datasetName]]$actual
    return(relevantAssignments[, calcForecastError(
      newDatasetPrev = newDatasetList[[datasetName]][newPrevIdx],
      newDatasetActual = newDatasetList[[datasetName]][newActIdx],
      oldDatasetPrev = oldDatasetList[[datasetName]][oldPrevIdx],
      oldDatasetActual = oldDatasetList[[datasetName]][oldActIdx],
      newDataAssignments = as.numeric(mget(newDataNewAssCols)),
      oldDataAssignments = as.numeric(mget(oldDataAssCols)),
      dissMatrix = dissMatrix, dissWeightConversion = dissWeightConversion),
      by = nonAssignmentColumns])
  })
  stopCluster(computingCluster)
  # Quick sanity check: for each clustering there should be a forecast error using
  # multiple forecast methods
  if (nrow(result) %% newDataRowCount != 0) {
    stop("Some cluster assignments are related to no dissimilarity matrix.")
  }
  result <- formatAndSortClusterSummaryTable(result)
  result[, Index := as.factor(Index)]
  result[, Reference := as.factor(Reference)]
  return(formatAndSortClusterSummaryTable(result))
}

# Summarizes forecasting errors on aggregated datasets with and without clustering:
# - takes two lists of list of paths pointing to the aggregated datasets (each
# element should represent a directory containing *.rds files, both lists should
# have the same element names and the same aggregations per datasets)
# - takes two lists of lists of two index vectors ("prev" and "actual") which
# determine which time series act as predecessors in the forecasting process and
# which time series act as exemplars for the target to forecast; the outer list
# should be named like the datasets; will be used for sub-setting the time series
# lists (input part, target part in forecasting)
# - takes two data.tables containing columns "Dataset", "Clustering", Dissimilarity",
# "Aggregation", "Level" and the clustering assignments (all remaining columns);
# "oldDataAssignmentsTable" should be the result of "analyzeAggregationImpact()",
# "newDataAssignmentsTable" should be the result of "assignDatasetsToClusters()"
# - takes a vector containing paths to dissimilarity matrices, describing the
# dissimilarity of the test time series to the training time series (returned
# by the assignment routine "assignDatasetsToClusters()")
# - creates a data.table with the same initial columns as the assignment tables
# and additional columns "Index", "Reference", "Value", containing the forecasting
# errors on the aggregated datasets using/not using the cluster assignments and
# dissimilarities for weighting
# - if "forecastAggregation" is provided, we still use dissimilarities of each
# aggregation type for weighting, but make the forecast itself always on the
# aggregation type "forecastAggregation" at each level (and fall back to forecast
# individual aggregations if the type is not available on a certain level); e.g.,
# you might want to make forecasts always with the mean series
createRelativeForecastErrorTable <- function(newDatasetBasePathList,
    oldDatasetBasePathList, newDataPrevActList, oldDataPrevActList,
    newDataAssignmentsTable, oldDataAssignmentsTable, dissMatrixPaths,
    dissWeightConversion = "minmax", forecastAggregation = NULL) {
  dissMatrixPaths <- list.files(dissMatrixPaths, pattern = "^DissMatrix_.*\\.rds$",
                                full.names = TRUE)
  # Change column types (should be same in both tables, avoid factors)
  newDataAssignmentsTable <- formatClusterSummaryTableForJoin(newDataAssignmentsTable)
  oldDataAssignmentsTable <- formatClusterSummaryTableForJoin(oldDataAssignmentsTable)
  # Create a joint assignmensTable because old-new pairs are used below
  # code the same as in createBaseForecastErrorTable()
  oldDataAssCols <- grep(pattern = "^A[0-9]+", x = colnames(oldDataAssignmentsTable),
                         value = TRUE)
  newDataCols <- colnames(newDataAssignmentsTable)
  newDataOldAssCols <- grep(pattern = "^A[0-9]+", x = newDataCols, value = TRUE)
  newDataNewAssCols <- gsub("^A", "A\\.", newDataOldAssCols)
  nonAssignmentColumns <- setdiff(newDataCols, newDataOldAssCols)
  # Rename so we can distinguish old and new assignments
  setnames(newDataAssignmentsTable, old = newDataOldAssCols, new = newDataNewAssCols)
  newDataRowCount <- nrow(newDataAssignmentsTable)
  newDataAssignmentsTable <- merge(newDataAssignmentsTable, oldDataAssignmentsTable,
                                   all.x = TRUE)
  if (nrow(newDataAssignmentsTable) != newDataRowCount) {
    stop("There is a Dataset-Clustering-Dissimilarity-Aggregation-Level combination
         in newDataAssignmentsTable which does not exist in oldDataAssignmentsTable.")
  }
  # Summarize datasets found under base paths, check further assumptions
  datasetOverviewTable <- checkAssignmentAssumptions(newDatasetBasePathList,
      oldDatasetBasePathList, newDataAssignmentsTable)
  # Forecasts using assignments, corresponding datasets and dissimilarity matrices
  computingCluster <- makeCluster(detectCores())
  registerDoParallel(computingCluster)
  result <- rbindlist(foreach(i = 1:nrow(datasetOverviewTable),
        .packages = c("data.table", "FastTSDistances"),
        .export = c("calcMAE", "calcForecastError")) %dopar% {
    currDataset <- as.character(datasetOverviewTable[i, Dataset])
    currAggregation <- as.character(datasetOverviewTable[i, Aggregation])
    currLevel <- datasetOverviewTable[i, Level]
    oldPath <- datasetOverviewTable[i, OldPath]
    newPath <- datasetOverviewTable[i, NewPath]
    if (!is.null(forecastAggregation)) {
      # Use dataset of a certain aggregation to make forecasts
      relevantOverviewTable <- datasetOverviewTable[Dataset == currDataset &
          Level == currLevel & Aggregation == forecastAggregation]
      if (nrow(relevantOverviewTable) == 1) {
        oldPath <- relevantOverviewTable[1, OldPath]
        newPath <- relevantOverviewTable[1, NewPath]
      } else {
        warning(paste0(
          "No datatset \"", currDataset, "\" with the level \"", currLevel, "\" ",
          "for the \"forecastAggregation\" parameter ",  "with the value \"",
          forecastAggregation, "\" found. Defaulting to the aggregation which ",
          "was used to create the clustering."))
      }
    }
    currOldDataset <- readRDS(oldPath)
    currNewDataset <- readRDS(newPath)
    # assumed naming convection (see analyzeAggregationImpactSingleDis()):
    # <<irrelevant>>_<<Dataset>>_<<Dissimilarity>>_<<Aggregation>>_<<Level>>.rds
    currDissMatrixPaths <- grep(paste0(currDataset, "_.+_",
        currAggregation, "_", currLevel, "\\.rds"), dissMatrixPaths, value = TRUE)
    if (length(currDissMatrixPaths) == 0) {
      stop(paste0("No dissimilarity matrix for ", currDataset, "-", currAggregation,
                  "-", currLevel))
    }
    return(rbindlist(lapply(seq_len(length(currDissMatrixPaths)), function(j) {
      # following code quite similar to createBaseForecastErrorTable()
      currDissimilarity <- strsplit(gsub("\\.rds", "", basename(currDissMatrixPaths[j])),
                                    "_")[[1]][3]
      relevantAssignments <- newDataAssignmentsTable[Dataset == currDataset &
          Dissimilarity == currDissimilarity & Aggregation == currAggregation &
          Level == currLevel]
      if (nrow(relevantAssignments) == 0) {
        stop(paste0("No suitable entry in newDataAssignmentsTable for dissimilarity ",
                    "matrix ", currDissMatrixPaths[j]))
      }
      currDissMatrix <- readRDS(currDissMatrixPaths[j])
      # Forecasts are made based on time series and their predecessors
      newPrevIdx <- newDataPrevActList[[currDataset]]$prev
      newActIdx <- newDataPrevActList[[currDataset]]$actual
      oldPrevIdx <- oldDataPrevActList[[currDataset]]$prev
      oldActIdx <- oldDataPrevActList[[currDataset]]$actual
      return(relevantAssignments[, calcForecastError(
        newDatasetPrev = currNewDataset[newPrevIdx],
        newDatasetActual = currNewDataset[newActIdx],
        oldDatasetPrev = currOldDataset[oldPrevIdx],
        oldDatasetActual = currOldDataset[oldActIdx],
        newDataAssignments = as.numeric(mget(newDataNewAssCols)),
        oldDataAssignments = as.numeric(mget(oldDataAssCols)),
        dissMatrix = currDissMatrix, dissWeightConversion = dissWeightConversion),
        by = nonAssignmentColumns])
    })))
  })
  stopCluster(computingCluster)
  # Quick sanity check: for each clustering there should be a forecast error using
  # multiple forecast methods
  if (nrow(result) %% newDataRowCount != 0) {
    stop("Some cluster assignments are related to no dissimilarity matrix
         and/or new/old dataset.")
  }
  result <- formatAndSortClusterSummaryTable(result)
  result[, Index := as.factor(Index)]
  result[, Reference := as.factor(Reference)]
  return(formatAndSortClusterSummaryTable(result))
}

# Takes a base/relative forecast error table and replaces the values for each
# forecasting strategy ("Reference" column) with the ratio of this value and
# the value of a "baseRef" reference for the same Dataset, Dissimilarity, Clustering,
# Aggregation, Level
transformToErrorRatioTable <- function(forecastErrorTable, baseRef = "naive") {
  baseRefTable <- forecastErrorTable[Reference == baseRef]
  baseRefTable[, Reference := NULL]
  setnames(baseRefTable, "Value", "BaseRefValue")
  result <- merge(forecastErrorTable, baseRefTable)
  result[, Value := Value / BaseRefValue]
  result[is.na(Value), Value := 1] # both techniques predict with error of 0
  result[, BaseRefValue := NULL]
  result <- formatAndSortClusterSummaryTable(result)
  return(result)
}

# Prints the number of unexpected NAs in a forecast error table.
checkForecastErrorNAs <- function(forecastErrorTable) {
  # no NAs expected
  cat(forecastErrorTable[, sum(is.na(Value))], "unexpected NAs for forecast error.\n")
  return(TRUE)
}

### Memory consumption ###

# Checks the files of all dissimilarity matrices on the "dissMatrixPath" and returns
# a table similar to other clustering evaluation tables, containing the sizes in the
# column "Value"; "inMemory" determines if we really load all matrices and check their
# memory consumption (slow) or if we take the file size on disk
summarizeDissMatrixMemoryConsumption <- function(dissMatrixPath, inMemory = FALSE) {
  dissMatrixNameParts <- strsplit(gsub("\\.rds", "", list.files(dissMatrixPath)), split = "_")
  if (inMemory) {# load and check object size
    dissMatrixSizes <- sapply(list.files(dissMatrixPath, full.names = TRUE),
                              function(filePath) object.size(readRDS(filePath)))
  } else {# on disk
    dissMatrixSizes <- file.size(list.files(dissMatrixPath, full.names = TRUE))
  }
  result <- rbindlist(lapply(seq_len(length(dissMatrixNameParts)), function(i) {
    list(Dataset = dissMatrixNameParts[[i]][2],
         Dissimilarity = dissMatrixNameParts[[i]][3],
         Aggregation = dissMatrixNameParts[[i]][4],
         Level = dissMatrixNameParts[[i]][5],
         Value = dissMatrixSizes[i])
  }))
  result[, Clustering := "test"] # "Clustering" col necessary for formatting
  result <- formatAndSortClusterSummaryTable(result)
  result[, Clustering := NULL]
  return(result)
}

# Checks the files of all time series lists on the "tsListPath" and returns
# a table similar to other clustering evaluation tables, containing the sizes in the
# column "Value"; "inMemory" determines if we really load lists and check their
# memory consumption (slow) or if we take the file size on disk
summarizeTimeSeriesMemoryConsumption <- function(tsListPath, inMemory = FALSE) {
  tsListNameParts <- strsplit(gsub("\\.rds", "", list.files(tsListPath)), split = "_")
  if (inMemory) {# load and check object size
    tsListSizes <- sapply(list.files(tsListPath, full.names = TRUE),
                              function(filePath) object.size(readRDS(filePath)))
  } else {# on disk
    tsListSizes <- file.size(list.files(tsListPath, full.names = TRUE))
  }
  result <- rbindlist(lapply(seq_len(length(tsListNameParts)), function(i) {
    list(Dataset = tsListNameParts[[i]][1],
         Aggregation = tsListNameParts[[i]][2],
         Level = tsListNameParts[[i]][3],
         Value = tsListSizes[i])
  }))
  result[, Clustering := "test"] # "Clustering" col necessary for formatting
  result[, Dissimilarity := "test"] # "Dissimilarity" col necessary for formatting
  result <- formatAndSortClusterSummaryTable(result)
  result[, c("Clustering", "Dissimilarity") := NULL]
  return(result)
}

### Plot routines ###

# Calculates breaks for a plot axis based on the limits (min, max). Creates
# 10-20 breaks.
calcEvaluationPlotBreaks <- function(limits) {
  powof10 <- floor(log10(limits[2] - limits[1]))
  multOfPowOf10 <- (limits[2] - limits[1]) / 10^powof10
  if (multOfPowOf10 <= 2) {
    breakDiff <- 10^(powof10 - 1)
  } else if (multOfPowOf10 <= 4) {
    breakDiff <- 2 * 10^(powof10 - 1)
  } else if (multOfPowOf10 <= 5) {
    breakDiff <- 2.5 * 10^(powof10 - 1)
  } else {
    breakDiff <- 5 * 10^(powof10 - 1)
  }
  return(seq(from = ceiling(limits[1] / breakDiff) * breakDiff,
             to = floor(limits[2] / breakDiff) * breakDiff, by = breakDiff))
}

# Plots the cluster validity indices from an "evaluationTable" (produced by
# "createClusteringEvaluationTable()"); different lines are created based on the
# "groupExpression" (which could be a column name or also a combination)
plotOneClusteringEvaluation <- function(evaluationTable, groupExpression = "Index",
                                        plotTitle = "") {
  evaluationTable <- evaluationTable[order(-Level)]
  # Convert Level to factor
  aggregationLevels <- evaluationTable[, sort(unique(Level), decreasing = TRUE)]
  evaluationTable[, Level := factor(Level, levels = aggregationLevels)]
  # External CVIs have a reference (e.g. previous clustering, base clustering), other
  # evaluation measures are only related to the clusterings independently
  hasReferenceColumn <- "Reference" %in% colnames(evaluationTable)
  # Find last level values (might be different in the groups), so we can put
  # group labels into plot (besides legend)
  lastLevelTable <- merge(evaluationTable,
                          evaluationTable[, .(Level = Level[.N]), by = groupExpression],
                          by = c(groupExpression, "Level"))
  if (hasReferenceColumn) {
    # labels preferrably to Reference == B, otherwise C, otherwise ...)
    lastLevelTable <- lastLevelTable[Reference == intersect(
      c("B", "C", "P", "M", "clus.weighted", "all.weighted", "1NN", "naive"),
      unique(Reference))[1]]
    plotLines <- lapply(evaluationTable[, unique(Reference)], function(refName) {
      geom_line(data = evaluationTable[Reference == refName],
                aes_string(x = "Level", y = "Value", color = groupExpression,
                    group = groupExpression, linetype = shQuote(refName)),
                na.rm = TRUE)
    })
  } else {
    plotLines <- list(
      geom_line(data = evaluationTable, aes(x = Level, y = Value,
        color = eval(parse(text = groupExpression)),
        group = eval(parse(text = groupExpression))))
    )
  }
  ggplot() +
    plotLines +
    scale_y_continuous(breaks = calcEvaluationPlotBreaks, minor_breaks = NULL) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.spacing = unit(0, "pt")) +
    scale_linetype_manual("Reference",
        values = c("B" = "solid", "C" = "dashed", "P" = "dashed", "M" = "twodash",
                   "clus" = "solid", "all" = "dashed", "1NN" = "twodash",
                   "naive" = "dotted"),
        labels = c("B" = "Base", "C" = "Current", "P" = "Previous", "M" = "Median",
                   "clus" = "clustered", "all" = "all", "1NN" = "1NN",
                   "naive" = "naive")) +
    scale_color_discrete(groupExpression) +
    geom_text(data = lastLevelTable, # additionally label lines at end
              mapping = aes(x = Level, y = Value, label = get(groupExpression),
                            color = get(groupExpression)), hjust = "inward", show.legend = FALSE) +
    ggtitle(label = plotTitle)
}

# Creates files showing how certain indices change with the aggregation level
# - "evaluationTable": should be created with one of the "create...Table()" methods,
# contains the raw data for the plots; apart from the attributes which are
# displayed in each plot, all remaining column-value combinations result in
# multiple plots
# - "comparisonAttribute": each plot will contain lines with different colors;
# you can choose what to group (e.g. different validity indices, dissimilarities etc.)
# - "plotPath": target directory for plots (created if not existing)
# - "plotName": name prefix for plot files
plotClusteringEvaluationsToFile <- function(evaluationTable, comparisonAttribute = "Index",
    plotPath = "./", plotName = "AggregationComparison", progressBar = TRUE) {
  if (!dir.exists(plotPath)) {
    dir.create(plotPath, recursive = TRUE)
  }
  # Flexible are all columns except "Level" (x-axis), "Value" (y-axis),
  # "Reference" (linetype) and "comparisonAttribute" (group in ggplot aes());
  # create plots for all value combinations of the "relevantColumns"
  relevantColumns <- setdiff(colnames(evaluationTable),
                             c("Level", "Value", "Reference", comparisonAttribute))
  relevantColumnValueCombis <- evaluationTable[, .GRP, by = relevantColumns][, -"GRP"]
  relevantColumnValueCombis <- relevantColumnValueCombis[, lapply(.SD, as.character)]
  plotCount <- nrow(relevantColumnValueCombis)
  computingCluster <- makeCluster(detectCores())
  registerDoParallel(computingCluster)
  foreach(i = 1:plotCount, .packages = c("data.table", "ggplot2", "tcltk"),
          .export = c("plotOneClusteringEvaluation", "calcEvaluationPlotBreaks")) %dopar% {
    if (progressBar) {
      if (!exists("progBar")) {
        progBar <- tkProgressBar(max = plotCount, title = "Plot progress")
      }
      setTkProgressBar(progBar, value = i, label = paste0("Creating plot ", i, "/", plotCount))
    }
    plotTitle <- paste(colnames(relevantColumnValueCombis), relevantColumnValueCombis[i],
                       sep = ": ", collapse = " - ")
    png(filename = paste0(plotPath, plotName, "_",
                          paste(relevantColumnValueCombis[i], collapse = "_"), ".png"),
        width = 960, height = 540)
    print(plotOneClusteringEvaluation(evaluationTable = merge(relevantColumnValueCombis[i],
          evaluationTable), groupExpression = comparisonAttribute, plotTitle = plotTitle))
    dev.off()
    return(TRUE)
  }
  stopCluster(computingCluster)
  return(TRUE)
}

# Plots all evaluation index values from an "evaluationTable" (produced by
# "createClusteringEvaluationTable()") against the aggregation level; if you
# provide a "groupExpression" (like Dataset, Aggregation, Clustering, Dissimilarity),
# it will be colored
plotClusteringEvaluationOverview <- function(evaluationTable, groupExpression = NA,
                                             plotTitle = "") {
  evaluationTable <- evaluationTable[!is.na(Value)]
  evaluationTable <- evaluationTable[order(-Level)]
  evaluationTable[, Level := as.factor(Level)]
  if (!is.na(groupExpression) && groupExpression %in% names(evaluationTable)) {
    ggplot() +
      geom_point(data = evaluationTable, mapping = aes_string(x = "Level",
          y = "Value", color = groupExpression), position = "jitter", size = 1) +
      ggtitle(plotTitle)
  } else {
    # geom_violin() is alternative to boxplot which shows density;
    # geom_quantile() draws regression lines over y quantiles at x values
    ggplot(data = evaluationTable, mapping = aes(x = Level, y = Value)) +
      geom_boxplot(color = "red", size = 2) +
      geom_point(position = "jitter", size = 1) +
      ggtitle(plotTitle)
  }
}
