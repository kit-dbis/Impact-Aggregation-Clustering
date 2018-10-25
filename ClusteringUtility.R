library(apcluster)
library(batchtools)
library(cluster)
library(dbscan)
library(data.table)
library(ggplot2)
library(FastTSDistances)
library(snow) # necessary for batchtools, but not imported there

source("GeneralPlotRoutinesUtility.R")

options(expressions = 100000) # "Error: evaluation nested too deeply" might occur otherwise

### Clustering algorithms with pre-configured parameters ###

# to work with do.call routines in the analysis, new clustering routines need to
# have following parameters; tsList, dissMatrix, tsPlotPath, clusPlotPath; further
# parameters are possible and should be listed here
DEFAULT_CLUSTERINGS <- list(
  AffProp = list(
    name = "AffProp", method = "clusterTSWithAffinityPropagation", params = list()
  ),
  DBSCAN = list(
    name = "DBSCAN", method = "clusterTSWithDBSCAN", params = list()
  ),
  Hier.avg = list(
    name = "Hier.avg", method = "clusterTSWithHierarchical",
    params = list(k = 2:10, linkageMethod = "average")
  ),
  Hier.com = list(
    name = "Hier.com", method = "clusterTSWithHierarchical",
    params = list(k = 2:10, linkageMethod = "complete")
  ),
  Hier.sin = list(
    name = "Hier.sin", method = "clusterTSWithHierarchical",
    params = list(k = 2:10, linkageMethod = "single")
  ),
  Hier.ward = list(
    name = "Hier.ward", method = "clusterTSWithHierarchical",
    params = list(k = 2:10, linkageMethod = "ward.D2")
  ),
  PAM = list(
    name = "PAM", method = "clusterTSWithPAM", params = list(k = 2:10)
  )
)

### Cluster plots ###

# Takes plots of time series from a "sourceDirectory" (which should have the "plotNames")
# and a vector of "clusterAssignments" (same length), creating subdirectories for each
# distinct assignment value and copies all corresponding time series plots
copyPlotsAccordingToClustering <- function(plotNames, clusterAssignments,
                                           sourceDirectory = "./", targetDirectory = "./") {
  stopifnot(length(plotNames) == length(clusterAssignments))
  sapply(unique(clusterAssignments), function(x) {
    dir.create(path = paste0(targetDirectory, "Cluster_", x),
               recursive = TRUE, showWarnings = FALSE)
  })
  for (i in 1:length(plotNames)) {
    file.copy(from = paste0(sourceDirectory, plotNames[i]),
              to = paste0(targetDirectory, "Cluster_", clusterAssignments[i], "/", plotNames[i]),
              overwrite = TRUE)
  }
}

# Creates a plot of all time series in a cluster, highlighting if (existing) the
# centroid and the standard deviatation at each item nr (only works for time
# series with equal length)
createClusterPlot <- function(tsMemberList, tsCentroid = NULL, plotTitle = "Cluster") {
  if (is.null(tsCentroid)) {
    return(
      ggplot() +
        lapply(tsMemberList, function(memberTS) {
          if (is.matrix(memberTS)) {#multivariate time series
            dataTable <- data.table(Item = seq_len(nrow(memberTS)), memberTS)
            dataTable <- melt(dataTable, id.vars = "Item", variable.name = "Attribute",
                              value.name = "Value")
            geom_dispatch(data = dataTable, aes(x = Item, y = Value, color = Attribute))
          } else {#univariate time series
            geom_dispatch(data = data.table(Item = seq_len(length(memberTS)),
                                        Value = memberTS),
                      mapping = aes(x = Item, y = Value), color = "black")
          }}) +
        ggtitle(plotTitle)
    )
  } else if (is.matrix(tsCentroid)) {# multivariate time series
    centroidTable <- data.table(Item = seq_len(nrow(tsCentroid)), tsCentroid)
    centroidTable <- melt(centroidTable, id.vars = "Item", variable.name = "Attribute",
                      value.name = "Value")
    if (length(tsMemberList) > 1 &&
        all(sapply(tsMemberList, nrow) == nrow(tsCentroid))) { # all same length
      if (nrow(tsCentroid) == 1) {
        sdMatrix <- matrix(apply(sapply(tsMemberList, identity), 1, PSDAA_fast, 1), nrow = 1)
      } else {
        sdMatrix <- sapply(1:ncol(tsCentroid), function(j) {
          apply(sapply(tsMemberList, function(ts) ts[, j]), 1, PSDAA_fast, 1)
        })
      }
      upperSDTable <- data.table(Item = seq_len(nrow(tsCentroid)), tsCentroid + sdMatrix)
      upperSDTable <- melt(upperSDTable, id.vars = "Item", variable.name = "Attribute",
                           value.name = "Value")
      lowerSDTable <- data.table(Item = seq_len(nrow(tsCentroid)), tsCentroid - sdMatrix)
      lowerSDTable <- melt(lowerSDTable, id.vars = "Item", variable.name = "Attribute",
                           value.name = "Value")
    } else {
      upperSDTable <- data.table(Item = integer(), Attribute = character(),
                                 Value = numeric())
      lowerSDTable <- upperSDTable
    }
    return(
      ggplot() +
        lapply(tsMemberList, function(memberTS) {
          dataTable <- data.table(Item = seq_len(nrow(memberTS)), memberTS)
          dataTable <- melt(dataTable, id.vars = "Item", variable.name = "Attribute",
                            value.name = "Value")
          geom_dispatch(data = dataTable, aes(x = Item, y = Value, color = Attribute),
                    alpha = .4)}) +
        geom_dispatch(data = centroidTable, aes(x = Item, y = Value, color = Attribute),
                  size = 1) +
        geom_dispatch(data = upperSDTable, aes(x = Item, y = Value, color = Attribute),
                  size = 1, linetype = "dashed") +
        geom_dispatch(data = lowerSDTable, aes(x = Item, y = Value, color = Attribute),
                  size = 1, linetype = "dashed") +
        ggtitle(plotTitle)
    )
  } else {# univariate time series
    if (length(tsMemberList) > 1 &&
        all(sapply(tsMemberList, length) == length(tsCentroid))) { # same length
      if (length(tsCentroid) == 1) {
        sdVec <- PSDAA_fast(as.numeric(tsMemberList), windowCount = 1)
      } else {
        sdVec <- apply(sapply(tsMemberList, identity), 1, PSDAA_fast, 1)
      }
      upperSDTable <- data.table(Item = seq_len(length(tsCentroid)),
                                 Value = tsCentroid + sdVec)
      lowerSDTable <- data.table(Item = seq_len(length(tsCentroid)),
                                 Value = tsCentroid - sdVec)
    } else {
      upperSDTable <- data.table(Item = integer(), Value = numeric())
      lowerSDTable <- upperSDTable
    }
    return(
      ggplot() +
        lapply(tsMemberList, function(memberTS) {
          geom_dispatch(data = data.table(Item = seq_len(length(memberTS)), Value = memberTS),
                        mapping = aes(x = Item, y = Value), color = "darkgrey")}) +
        geom_dispatch(data = data.table(Item = seq_len(length(tsCentroid)), Value = tsCentroid),
                  mapping = aes(x = Item, y = Value), color = "red") +
        geom_dispatch(data = lowerSDTable, mapping = aes(x = Item, y = Value), color = "black") +
        geom_dispatch(data = upperSDTable, mapping = aes(x = Item, y = Value), color = "black") +
        ggtitle(plotTitle)
    )
  }
}

# Creates plots for a clustering results and saves them on the hard drive. If
# "clusPlotPath" is set one plot will be created for each cluster depicting
# all time series in it and also the centroid (if existing). If "tsPlotPath" is
# set, the time series plots located there will also be copied into sub-directories
# of the "clusPlotPath" based on the cluster assignments.
# If you pass centroid ids, make sure that there are as many as distinct clusters
# and the centroid ids are in the same order as the sorted cluster labels (which
# should usually be fullfilled if you have assignments in 1:k and a k-length vector
# of centroids)
saveClusterPlots <- function(tsList, assignments, centroidIdx = NULL, tsPlotPath = NULL,
                             clusPlotPath = NULL) {
  if (!is.null(clusPlotPath)) {
    dir.create(clusPlotPath, showWarnings = FALSE, recursive = TRUE)
    clusterIds <- sort(unique(assignments))
    if (!is.null(centroidIdx) && length(clusterIds) != length(centroidIdx)) {
      stop("There should be as many centroids as clusters.")
    }
    for (i in seq_len(length(clusterIds))) {
      if (is.null(centroidIdx)) {
        tsCentroid <- NULL
      } else {
        tsCentroid <- tsList[[centroidIdx[i]]]
      }
      png(filename = paste0(clusPlotPath, "Cluster_", i, ".png"),
          width = 960, height = 540)
      print(createClusterPlot(tsMemberList = tsList[assignments == clusterIds[i]],
                              tsCentroid = tsCentroid, plotTitle = paste0("Cluster", i)))
      dev.off()
    }
    if (!is.null(tsPlotPath)) {
      copyPlotsAccordingToClustering(plotNames = paste0(names(tsList), ".png"),
          clusterAssignments = assignments, sourceDirectory = tsPlotPath,
          targetDirectory = clusPlotPath)
    }
  }
  return(TRUE)
}

### Clustering + validation routines ###

# Calculates the Silhouette coefficient; returns NA if there is only one
# cluster or there are only singleton clusters
avgSilhouette <- function(clusterAssignments, dissMatrix) {
  # silhouette() cannot handle Inf values (as produced by CID dissimilarities);
  # replace with very high values instead (which muss be small enough that sums
  # are smaller than Inf, too -> divide by max cluster size, account for numerical
  # errors by a factor of two)
  replacementValue <- .Machine$double.xmax / (2 * max(tabulate(clusterAssignments)))
  dissMatrix[dissMatrix > replacementValue] <- replacementValue
  silObject <- silhouette(clusterAssignments, dist = dissMatrix)
  if (length(silObject) == 1 && is.na(silObject)) {
    return(NA)
  } else {
    return(mean(silObject[, "sil_width"]))
  }
}

# Computes the connectivity measure of *Handl 2005*, but normalizes it to [0,1]
# (connectivity can be at most (1 + 1/2 + 2/3 + 1/k)*n for kNN)
# and also inverts it such that high values are good; k = number of neighbors
normConnectivity <- function(clusterAssignments, dissMatrix, k = 10) {
  return(1 - connectivity.diss.mx(diss.mx = dissMatrix, clust = clusterAssignments,
      neighbour.num = k) / (length(clusterAssignments) * sum(1/(1:k))))
}

# Use the affinity propagation algorithm [package: "apcluster"] to cluster a list
# of time series; parameters:
# - distance matrix has to be pre-computed and passed, list of time series
# is only used to plot centroids
# - an input path to the plotted time series and an output path for the clustering
# result, so the plots will be copied according to the clustering and centroid
# plots will be created
# - no parameters of clustering algorithm itself, reasonable defaults chosen
clusterTSWithAffinityPropagation <- function(tsList, dissMatrix, tsPlotPath = NULL,
    clusPlotPath = NULL) {
  cat("Clustering with Affinity Propagation.\n")
  # Dissimilarity to similarity conversion: just put minus in front (apcluster can handle this)
  # Input preference: none (p=NA), median of values in s (q=NA / q = .5) [default/paper]
  # Iterations (absolute maximum / earlier stop if no change): 1000/100 [default]
  # Dampening factor "lam": 0.9 [default, original paper uses weaker 0.5]
  # Add small amount of noise for stability [default]
  # Set random seed (only for noise)
  apropResult <- apcluster(s = -dissMatrix, p = NA, q = .5, maxits = 1000,
      convits = 100, lam = .9, nonoise = FALSE, seed = 25)
  apropAssignments <- labels(apropResult, type = "enum")
  saveClusterPlots(tsList = tsList, assignments = apropAssignments,
      centroidIdx = apropResult@exemplars, tsPlotPath = tsPlotPath,
      clusPlotPath = clusPlotPath)
  return(apropAssignments)
}

# Use the DBSCAN algorithm to cluster a list of time series; parameters:
# - distance matrix has to be pre-computed and passed, list of time series
# is only used to plot centroids
# - an input path to the plotted time series and an output path for the clustering
# result, so the plots will be copied according to the clustering and centroid
# plots will be created
# - minPts set to 1 so that there are no noise points (which would cause problems
# in all evaluation routines, the is no established approach how to deal with noise
# in validation, *Moulavi 2014*: "Density-Based Clustering Validation")
# - epsilon determined as mean 1NN distance (very simple heuristic); many standard
# validity indices favour results with only one big cluster and a few or even no
# singletons, so there is no good way for parameter choice by comparison as we do
# for PAM and hclust (approach of *Moulavi 2014* quite complicated)
clusterTSWithDBSCAN <- function(tsList, dissMatrix, tsPlotPath = NULL,
    clusPlotPath = NULL) {
  cat("Clustering with DBSCAN.\n")
  dissObject <- as.dist(dissMatrix)
  diag(dissMatrix) <- Inf
  epsilon <- mean(apply(dissMatrix, 1, min)) # mean 1NN distance
  dbscanResult <- dbscan(dissObject, eps = epsilon, minPts = 1)$cluster
  saveClusterPlots(tsList = tsList, assignments = dbscanResult,
                   tsPlotPath = tsPlotPath, clusPlotPath = clusPlotPath)
  return(dbscanResult)
}

# Use the hclust (agglomerative hierarchical) algorithm to cluster a list of time
# series; parameters:
# - distance matrix has to be pre-computed and passed, list of time series
# is only used to plot centroids
# - an input path to the plotted time series and an output path for the clustering
# result, so the plots will be copied according to the clustering and centroid
# plots will be created
# - number of clusters k, can also be a vector; the silhouette coefficient
# (average silhouette of all elements) will be used to determine the optimal
# clustering in that case
# - linkage method for agglomerative merges (passed to "hclust")
clusterTSWithHierarchical <- function(tsList, dissMatrix, tsPlotPath = NULL,
    clusPlotPath = NULL, k = 10, linkageMethod = "average") {
  cat("Clustering with hierarchical and", linkageMethod, "linkage.\n")
  dissMatrix[is.infinite(dissMatrix)] <- .Machine$double.xmax # hclust cannot handle Inf
  hclustResult <- hclust(d = as.dist(dissMatrix), method = linkageMethod)
  hclustResult <- cutree(hclustResult, k = k)
  if (length(k) > 1) { # each column of hclustResult represents on k
    silhouettes <- apply(hclustResult, 2, function(x)
      avgSilhouette(clusterAssignments = x, dissMatrix = dissMatrix))
    bestK <- which.max(silhouettes)
    if (length(bestK) == 0) { # all silhouettes NA (NAs/Inf in dissMatrix)
      bestK <- k[1]
    }
    hclustResult <- hclustResult[, bestK]
  }
  saveClusterPlots(tsList = tsList, assignments = hclustResult,
                   tsPlotPath = tsPlotPath, clusPlotPath = clusPlotPath)
  return(hclustResult)
}

# Use the PAM algorithm to cluster a list of time series; parameters:
# - distance matrix has to be pre-computed and passed, list of time series
# is only used to plot centroids
# - an input path to the plotted time series and an output path for the clustering
# result, so the plots will be copied according to the clustering and centroid
# plots will be created
# - number of clusters k, can also be a vector; the Silhouette coefficient
# (average Silhouette of all elements) will be used to determine the optimal
# clustering in that case
clusterTSWithPAM <- function(tsList, dissMatrix, tsPlotPath = NULL,
                             clusPlotPath = NULL, k = 10) {
  cat("Clustering with PAM.\n")
  if (length(k) == 1) {
    pamResult <- pam(x = dissMatrix, k = k, diss = TRUE)
  } else {
    pamResults <- lapply(k, function(numClusters) {
      return(pam(x = dissMatrix, k = numClusters, diss = TRUE))
    })
    pamSilhouettes <- sapply(pamResults, function(x) x$silinfo$avg.width)
    pamResult <- pamResults[[which.max(pamSilhouettes)]]
  }
  saveClusterPlots(tsList = tsList, assignments = pamResult$clustering,
                   centroidIdx = pamResult$id.med, tsPlotPath = tsPlotPath,
                   clusPlotPath = clusPlotPath)
  return(pamResult$clustering)
}

### Impact analysis ###

# Runs multiple clustering routines (stored in "clusList") on a dissimilarity
# matrix which will be pre-computed with the dissimilarity represented by "disEntry".
# Aggregation information is only necessary for creating a proper results list
# (which describes the general variables of the experiment and the cluster assignments).
# If "tsPlotPath" and "clusPlotBasePath" are supplied, time series plots are copied
# according to the clustering result from the "tsPlotPath" to folders in the
# "clusPlotBasePath" (assumes that plots have the same names as the elements in
# "tsList"). If "dissMatrixPath" is supplied, then the dissimilarity matrix
# will be saved (might be necessary for some evaluation routines),
analyzeAggregationImpactSingleDis <- function(tsList, datasetName, disEntry,
    clusList, aggName, aggLevel, dissMatrixPath = NULL, tsPlotPath = NULL,
    clusPlotBasePath = NULL, cpuCores = 8) {
  cat("Clustering dataset \"", datasetName, "\" aggregated as \"", aggName,
      "\" on level ", aggLevel, " with dissimilarity \"", disEntry$name, "\".\n",
      sep = "")
  clusteringResultsForDis <- list()
  # Dissimilarity matrix first (can be used for multiple clustering algorithms)
  startTime <- as.numeric(Sys.time())
  if (disEntry$takesList) {
    dissMatrix <- do.call(disEntry$method, args = c(list(tsList), disEntry$params))
  } else {
    dissMatrix <- sapply(seq_len(length(tsList)), function(j) { # cols
      if (j == 1) {
        return(c(0, rep(NA, times = length(tsList) - 1)))
      } else {
        return(c(sapply(1:(j - 1), function(i) { # rows
          return(do.call(what = disEntry$method, args = c(list(tsList[[i]],
              tsList[[j]]), disEntry$params)))
        }), 0, rep(NA, times = length(tsList) - j)))
      }
    })
    # Use symmetry of dissimilarities
    dissMatrix[lower.tri(dissMatrix, diag = FALSE)] <-
      t(dissMatrix)[lower.tri(dissMatrix, diag = FALSE)]
  }
  endTime <- as.numeric(Sys.time())
  dissComputationTime <- endTime - startTime
  if (!is.null(dissMatrixPath)) {
    # Be aware that naming scheme is used in evaluation routines, so only change
    # with care (or better not at all)
    saveRDS(dissMatrix, paste0(dissMatrixPath, "DissMatrix_", datasetName, "_",
                               disEntry$name, "_", aggName, "_", aggLevel, ".rds"))
  }
  # Clustering algorithms
  clusteringResultsForDis <- lapply(clusList, function(clusEntry) {
    if (is.null(clusPlotBasePath)) {
      clusPlotPath <- NULL
    } else {
      clusPlotPath <- paste0(clusPlotBasePath, disEntry$name, "_", clusEntry$name,
          "/", aggName, "_", aggLevel, "_Results/")
    }
    oneClusteringResult <- list(Dataset = datasetName, Dissimilarity = disEntry$name,
        Clustering = clusEntry$name, Aggregation = aggName, Level = aggLevel)
    startTime <- as.numeric(Sys.time())
    clusteringResult <- do.call(clusEntry$method,
        args = c(tsList = list(tsList), dissMatrix = list(dissMatrix),
                 tsPlotPath = tsPlotPath, clusPlotPath = clusPlotPath, clusEntry$params))
    endTime <- as.numeric(Sys.time())
    clusComputationTime <- endTime - startTime
    oneClusteringResult[paste0("A", 1:length(tsList))] <- clusteringResult
    oneClusteringResult[["DissimilarityTime"]] <- dissComputationTime
    oneClusteringResult[["ClusteringTime"]] <- clusComputationTime
    return(oneClusteringResult)
  })
  names(clusteringResultsForDis) <- paste0("Assignment_", datasetName, "_",
      disEntry$name, "_", sapply(clusList, function(clusEntry) clusEntry$name), 
      "_", aggName, "_", aggLevel)
  return(clusteringResultsForDis)
}

# Stops if assumptions of the aggregation-impact-on-clustering routines are violated,
# It is important that list elements are named, because the method selects datasets,
# dissimilarities etc. from lists according to names in the experiments' parameter
# tables.
# Common sub-routine for analysis based on raw as well as pre-aggregated data.
checkAnalysisAssumptions <- function(datasetList, disList, dissMatrixPathList,
    tsPlotPathList, clusPlotBasePathList) {
  if (is.null(names(datasetList)) ||
      any(sapply(names(datasetList), is.null)) ||
      any(sapply(names(datasetList), is.na)) ||
      any(names(datasetList) == "")) {
    stop("All your datasets should be named.")
  }
  if (is.null(names(disList)) ||
      any(sapply(names(disList), is.null)) ||
      any(sapply(names(disList), is.na)) ||
      any(names(disList) == "")) {
    stop("All your dissimilarities should be named.")
  }
  if (length(dissMatrixPathList) > 0 &&
      (is.null(names(dissMatrixPathList)) ||
       any(sapply(names(dissMatrixPathList), is.null)) ||
       any(sapply(names(dissMatrixPathList), is.na)) ||
       any(names(dissMatrixPathList) == "") ||
       length(setdiff(names(dissMatrixPathList), names(datasetList))) > 0)) {
    stop("All your paths for dissimilarity matrices should be named like the 
         corresponding dataset in your dataset list.")
  }
  if (length(tsPlotPathList) > 0 &&
      (is.null(names(tsPlotPathList)) ||
       any(sapply(names(tsPlotPathList), is.null)) ||
       any(sapply(names(tsPlotPathList), is.na)) ||
       any(names(tsPlotPathList) == "") ||
       length(setdiff(names(tsPlotPathList), names(datasetList))) > 0)) {
    stop("All your paths to time series plots should be named like the 
         corresponding dataset in your dataset list.")
  }
  if (length(clusPlotBasePathList) > 0 &&
      (is.null(names(clusPlotBasePathList)) ||
       any(sapply(names(clusPlotBasePathList), is.null)) ||
       any(sapply(names(clusPlotBasePathList), is.na)) ||
       any(names(clusPlotBasePathList) == "") ||
       length(setdiff(names(clusPlotBasePathList), names(datasetList))) > 0)) {
    stop("All your paths to clustering plots should be named like the 
         corresponding dataset in dataset list.")
  }
}

# Stops if assumptions of the aggregation-impact-on-clustering routine are violated,
# It is important that list elements are named, because the method selects datasets,
# dissimilarities etc. from lists according to names in the experiments' parameter
# tables.
# Routine for raw data which are aggregated during the analysis.
checkAnalysisAssumptions.raw <- function(tsListList, disList, aggList,
    dissMatrixPathList, tsPlotPathList, clusPlotBasePathList) {
  checkAnalysisAssumptions(datasetList = tsListList, disList = disList,
      dissMatrixPathList = dissMatrixPathList, tsPlotPathList = tsPlotPathList,
      clusPlotBasePathList = clusPlotBasePathList)
  if (is.null(names(aggList)) ||
      any(sapply(names(aggList), is.null)) ||
      any(sapply(names(aggList), is.na)) ||
      any(names(aggList)  == "")) {
    stop("All your aggregations should be named.")
  }
}

# Computes the cross-product of datasets, aggregations, aggregation levels,
# dissimilarities and clusterings using the "batchtools" package.
# Returns a data.table containing experiment parameters in the first columns
# and cluster assignments in the last columns.
# Routine is based on raw data, so aggregations are computed before dissimilarity
# computation and clustering start.
analyzeAggregationImpact.raw <- function(tsListList, disList, clusList, aggList,
    dissMatrixPathList = list(), tsPlotPathList = list(), clusPlotBasePathList = list(),
    expRegistry = NULL, cpuCores = 8) {
  checkAnalysisAssumptions.raw(tsListList = tsListList, disList = disList,
        aggList = aggList, dissMatrixPathList = dissMatrixPathList,
        tsPlotPathList = tsPlotPathList, clusPlotBasePathList = clusPlotBasePathList)
  if (is.null(expRegistry)) {
    expRegistry <- makeExperimentRegistry(file.dir = NA, seed = 1)
  }
  expRegistry$cluster.functions <- makeClusterFunctionsSocket(ncpus = cpuCores)
  expRegistry$source <- c("ClusteringUtility.R", "DissimilarityUtility.R")
  # Problems: different datasets combined with different aggregation and levels
  problemDescriptionList <- list()
  for (i in 1:length(tsListList)) {
    datasetName <- names(tsListList)[i]
    datasetParameterList <- list()
    for (j in 1:length(aggList)) {
      aggregationName <- names(aggList)[j]
      aggLevels <- do.call(
        what = aggList[[j]]$levelMethod,
        args = c(tsList = list(tsListList[[i]]), aggList[[j]]$levelParams)
      )
      datasetParameterList <- c(datasetParameterList, list(expand.grid(
        aggName = aggregationName, aggLevel = aggLevels
      )))
    }
    datasetParameterTable <- data.table(rbindlist(datasetParameterList))
    setorder(datasetParameterTable, -aggLevel)
    if (!is.null(dissMatrixPathList[[datasetName]])) {
      datasetParameterTable[, dissMatrixPath := dissMatrixPathList[[datasetName]]]
      dir.create(dissMatrixPathList[[datasetName]], showWarnings = FALSE, recursive = TRUE)
    }
    if (!is.null(tsPlotPathList[[datasetName]])) {
      datasetParameterTable[, tsPlotPath := tsPlotPathList[[datasetName]]]
    }
    if (!is.null(clusPlotBasePathList[[datasetName]])) {
      datasetParameterTable[, clusPlotBasePath := clusPlotBasePathList[[datasetName]]]
    }
    problemDescriptionList[[datasetName]] <- datasetParameterTable
    addProblem(
      name = datasetName,
      data = list(tsList = tsListList[[i]], disList = disList, clusList = clusList,
                  aggList = aggList),
      fun = function(data, job, ...) { # returns problem instance
        # Combine dataset + aggregation + level to get experiment datasets
        return(do.call(
          what = data$aggList[[job$pars$prob.pars$aggName]]$aggMethod,
          args = c(tsList = list(data$tsList),
                   aggLevel = job$pars$prob.pars$aggLevel,
                   data$aggList[[job$pars$prob.pars$aggName]]$aggParams)
        ))
      },
      reg = expRegistry
    )
  }
  # Algorithm, experiment, run, collect results
  return(analyzeAggregationImpact.algorithm(disList, problemDescriptionList, expRegistry))
}

# Computes the cross-product of datasets, aggregations, aggregation levels,
# dissimilarities and clusterings using the "batchtools" package.
# Returns a data.table containing experiment parameters in the first columns
# and cluster assignments in the last columns.
# This sub-routine expects a batchtools registry which already has the relevant
# "problems" (dataset and instances) added. The focus here is on building the
# dissimilarity-dependend "algorithms", add "experiments", run the jobs and
# collect the results.
analyzeAggregationImpact.algorithm <- function(disList, problemDescriptionList, expRegistry) {
  # Algorithms: different dissimilarities (same clusterings for each)
  addAlgorithm(name = "TSClustering", fun = function(data, job, instance, ...) {
    return(analyzeAggregationImpactSingleDis(
      tsList = instance,
      datasetName = job$prob.name,
      disEntry = data$disList[[job$pars$algo.pars$disName]],
      clusList = data$clusList,
      aggName = job$pars$prob.pars$aggName,
      aggLevel = job$pars$prob.pars$aggLevel,
      dissMatrixPath = job$pars$prob.pars$dissMatrixPath,
      tsPlotPath = job$pars$prob.pars$tsPlotPath,
      clusPlotBasePath = job$pars$prob.pars$clusPlotBasePath,
      cpuCores = 1
    )) # distance computations not parallel as jobs already are
  }, reg = expRegistry)
  # Experiment: Combine problems + algorithms
  addExperiments(
    prob.designs = problemDescriptionList,
    algo.designs = list(TSClustering = data.table(disName = names(disList))),
    reg = expRegistry
  )
  submitJobs(reg = expRegistry)
  waitForJobs(reg = expRegistry)
  # Results can be retrieved as list with one element per job, but each of these
  # elements summarizes multiple cluster assignments (as lists), therefore "unlist"
  return(rbindlist(unlist(reduceResultsList(reg = expRegistry), recursive = FALSE)))
}

# Stops if assumptions of the aggregation-impact-on-clustering routine are violated,
# It is important that list elements are named, because the method selects datasets,
# dissimilarities etc. from lists according to names in the experiments' parameter
# tables.
# Routine for analysis with pre-aggregated data from the hard-drive.
checkAnalysisAssumptions.preagg <- function(aggBasePathList, disList,
    dissMatrixPathList, tsPlotPathList, clusPlotBasePathList) {
  checkAnalysisAssumptions(datasetList = aggBasePathList, disList = disList,
      dissMatrixPathList = dissMatrixPathList, tsPlotPathList = tsPlotPathList,
      clusPlotBasePathList = clusPlotBasePathList)
  if (any(!dir.exists(unlist(aggBasePathList)))) {
    stop("All your paths to aggregated data files should exist.")
  }
}

# Computes the cross-product of datasets, aggregations, aggregation levels,
# dissimilarities and clusterings using the "batchtools" package.
# Returns a data.table containing experiment parameters in the first columns
# and cluster assignments in the last columns.
# Routine is based on pre-aggregated data which are retrieved from the hard drive
# before dissimilarity computation and clustering start.
analyzeAggregationImpact.preagg <- function(aggBasePathList, disList, clusList,
    dissMatrixPathList = list(), tsPlotPathList = list(), clusPlotBasePathList = list(),
    expRegistry = NULL, cpuCores = 8) {
  checkAnalysisAssumptions.preagg(aggBasePathList = aggBasePathList,
      disList = disList, dissMatrixPathList = dissMatrixPathList,
      tsPlotPathList = tsPlotPathList, clusPlotBasePathList = clusPlotBasePathList)
  if (is.null(expRegistry)) {
    expRegistry <- makeExperimentRegistry(file.dir = NA, seed = 1)
  }
  expRegistry$cluster.functions <- makeClusterFunctionsSocket(ncpus = cpuCores)
  expRegistry$source <- c("ClusteringUtility.R", "DissimilarityUtility.R")
  # Problems: different datasets combined with different aggregation and levels;
  # loaded from hard drive
  problemDescriptionList <- list()
  for (i in 1:length(aggBasePathList)) {
    datasetName <- names(aggBasePathList)[i]
    datasetBasePath <- aggBasePathList[[datasetName]]
    # Find aggregated data files (naming pattern -> AggregationUtility.R)
    datasetParameterTable <- data.table(rbindlist(lapply(
      list.files(path = datasetBasePath,
                 pattern = paste0(datasetName, "_.*_.*\\.rds")),
      function(dataFileName) {
        return(list(aggName = strsplit(dataFileName, split = "_")[[1]][[2]],
                    aggLevel = as.numeric(gsub(".rds", "", strsplit(dataFileName, split = "_")[[1]][[3]])),
                    aggFilePath = paste0(datasetBasePath, dataFileName)))
    })))
    setorder(datasetParameterTable, -aggLevel)
    if (!is.null(dissMatrixPathList[[datasetName]])) {
      datasetParameterTable[, dissMatrixPath := dissMatrixPathList[[datasetName]]]
      dir.create(dissMatrixPathList[[datasetName]], showWarnings = FALSE, recursive = TRUE)
    }
    if (!is.null(tsPlotPathList[[datasetName]])) {
      datasetParameterTable[, tsPlotPath := tsPlotPathList[[datasetName]]]
    }
    if (!is.null(clusPlotBasePathList[[datasetName]])) {
      datasetParameterTable[, clusPlotBasePath := clusPlotBasePathList[[datasetName]]]
    }
    problemDescriptionList[[datasetName]] <- datasetParameterTable
    addProblem(
      name = datasetName,
      data = list(disList = disList, clusList = clusList),
      fun = function(data, job, ...) { # returns problem instance
        # Combine dataset + aggregation + level to get experiment datasets
        return(readRDS(file = job$pars$prob.pars$aggFilePath))
      },
      reg = expRegistry
    )
  }
  # Algorithm, experiment, run, collect results
  return(analyzeAggregationImpact.algorithm(disList, problemDescriptionList, expRegistry))
}

### Cluster assignment ###

# Takes a list of directories (absolute, relative) and searches for .rds file
# representing aggregated datasets (which follow a certain naming scheme),
# returning a data.table with the columns "Dataset", "Aggregation", "Level"
summarizeAggregatedDatasets <- function(datasetBasePathList) {
  result <- rbindlist(lapply(datasetBasePathList, function(datsetBasePath) {
    # naming pattern -> AggregationUtility.R
    datasetAggFiles <- list.files(path = datsetBasePath, pattern = ".*_.*_.*\\.rds$")
    fileNameSplits <- strsplit(datasetAggFiles, split = "_")
    return(data.table(
      Dataset = as.factor(sapply(fileNameSplits, function(x) x[[1]])),
      Aggregation = as.factor(sapply(fileNameSplits, function(x) x[[2]])),
      Level = as.numeric(gsub(".rds", "", sapply(fileNameSplits, function(x) x[[3]]))),
      Path = paste0(datsetBasePath, datasetAggFiles)
    ))
  }))
  setorder(result, -Level)
  return(result)
}

# Makes sure that the base paths lists are named and have the same names. Guarantees
# that for each new dataset, there is an entry in the assignments table.
# Stops if an assumption is violated.
# Returns a data.table summarizing the detected old and new datasets (based on
# file names) together with their paths.
checkAssignmentAssumptions <- function(newDataBasePathList, oldDataBasePathList,
                                       clusterAssignmentsTable, disList = NULL) {
  if (!is.null(disList) &&
      (is.null(names(disList)) ||
       any(sapply(names(disList), is.null)) ||
       any(sapply(names(disList), is.na)) ||
       any(names(disList) == ""))) {
    stop("All your dissimilarities should be named.")
  }
  if (length(oldDataBasePathList) > 0 &&
      (is.null(names(oldDataBasePathList)) ||
       any(sapply(names(oldDataBasePathList), is.null)) ||
       any(sapply(names(oldDataBasePathList), is.na)) ||
       any(names(oldDataBasePathList) == ""))) {
    stop("All your base paths to the old data should be named.")
  }
  if (length(newDataBasePathList) > 0 &&
      (is.null(names(newDataBasePathList)) ||
       any(sapply(names(newDataBasePathList), is.null)) ||
       any(sapply(names(newDataBasePathList), is.na)) ||
       any(names(newDataBasePathList) == "") ||
       length(setdiff(names(newDataBasePathList), names(oldDataBasePathList))) > 0)) {
    stop("All your base paths to the new data should be named like the
         corresponding base paths to the old data.")
  }
  assignmentsMatrix <- as.matrix(clusterAssignmentsTable[, mget(grep("^A[0-9]+",
      colnames(clusterAssignmentsTable), value = TRUE))])
  if (!is.integer(assignmentsMatrix) ||
      any(assignmentsMatrix < 1)) {
    stop("The assignments should be positive integer values.")
  }
  aggNewDatasets <- summarizeAggregatedDatasets(newDataBasePathList)
  setnames(aggNewDatasets, old = "Path", new = "NewPath")
  # Join with assignment table to check if all Dataset-Aggregation-Level
  # combinations are covered there
  assignmentDatasets <- clusterAssignmentsTable[, 1, by = .(Dataset, Aggregation, Level)]
  aggDatasetsOverviewTable <- merge(x = aggNewDatasets, y = assignmentDatasets, all.x = TRUE)
  if (aggNewDatasets[, .N] != aggDatasetsOverviewTable[, .N]) {
    stop("There is one file under a newDataBasePath whose corresponding combination
         of Dataset, Aggregation and Level could not be found in the table of cluster
         assignments.")
  }
  # Join with corresponding table for old [clustered] datasets to check if all
  # Dataset-Aggregation-Level combinations are covered there
  aggOldDatasets <- summarizeAggregatedDatasets(oldDataBasePathList)
  setnames(aggOldDatasets, old = "Path", new = "OldPath")
  aggDatasetsOverviewTable <- merge(x = aggNewDatasets, y = aggOldDatasets, all.x = TRUE)
  if (aggNewDatasets[, .N] != aggDatasetsOverviewTable[, .N]) {
    stop("There is one file under a newDataBasePath whose corresponding dataset
         file in oldDataBasePathList could not be found.")
  }
  return(aggDatasetsOverviewTable)
}

# For each row in the "dissMatrix", the "numNeighbors" nearest neighbors are
# determined and the majority class from "oldAssignments" is returned (assuming
# that assignments are 1:k vectors with k = number of clusters); ties are
# resolved by taking the lower cluster number
determineAssignments <- function(dissMatrix, oldAssignments, numNeighbors) {
  if (numNeighbors == 1) {
    return(apply(dissMatrix, 1, function(dissRow) {
      oldAssignments[which.min(dissRow)]
    }))
  } else {
    return(apply(dissMatrix, 1, function(dissRow) {
      relevantNeighbors <- which(rank(dissRow) <= numNeighbors)
      # tabulate does only work for positive integers!
      which.max(tabulate(oldAssignments[relevantNeighbors]))
    }))
  }
}

# Takes two datasets, a clustered one and another one to be assigned, which are based
# on the same original dataset, aggregation type and level. Computes the distances
# from the new to the old time series and assigns cluster labels based on the
# "numNeighbors" nearest neighbors. Returns the assignments and the dissimilarity
# matrix in a list (as "analyzeAggregationImpactSingleDis()").
assignOneDatasetToClusters <- function(oldDataset, newDataset, oldAssignmentsList,
    numNeighbors, datasetName, disEntry, aggName, aggLevel, dissMatrixPath = NULL) {
  assignmentResultsForDis <- list()
  # Pre-compute dissimilarities (relevant for multiple clusterings)
  startTime <- as.numeric(Sys.time())
  if (disEntry$takesList) {
    dissMatrix <- do.call(disEntry$method, args = c(list(list(newDataset, oldDataset)),
                                                    disEntry$params))
  } else {
    dissMatrix <- matrix(nrow = length(newDataset), ncol = length(oldDataset))
    for (i in seq_len(length(newDataset))) {
      for (j in seq_len(length(oldDataset))) {
        dissMatrix[i,j] <- do.call(disEntry$method,
            args = c(list(newDataset[[i]], oldDataset[[j]]), disEntry$params))
      }
    }
  }
  endTime <- as.numeric(Sys.time())
  dissComputationTime <- endTime - startTime
  # Be aware that naming scheme of clustering results list is used in evaluation
  # routines, so only change with care (or better not at all)
  if (!is.null(dissMatrixPath)) {
    saveRDS(dissMatrix, paste0(dissMatrixPath, "DissMatrix_", datasetName, "_",
                               disEntry$name, "_", aggName, "_", aggLevel, ".rds"))
  }
  # Assign (determine best cluster) for each clustering result (element of oldAssignmentsList)
  assignmentResultsForDis <- lapply(seq_len(length(oldAssignmentsList)), function(i) {
    oneNewAssignment <- list(Dataset = datasetName, Clustering = names(oldAssignmentsList)[i],
      Dissimilarity = disEntry$name, Aggregation = aggName, Level = aggLevel)
    startTime <- as.numeric(Sys.time())
    assignmentResult <- determineAssignments(dissMatrix, oldAssignmentsList[[i]], numNeighbors)
    endTime <- as.numeric(Sys.time())
    assignmentComputationTime <- endTime - startTime
    oneNewAssignment[paste0("A", 1:length(newDataset))] <- assignmentResult
    oneNewAssignment[["DissimilarityTime"]] <- dissComputationTime
    oneNewAssignment[["AssignmentTime"]] <- assignmentComputationTime
    return(oneNewAssignment)
  })
  names(assignmentResultsForDis) <- paste0("Assignment_", datasetName, "_",
      names(oldAssignmentsList), "_", disEntry$name, "_", aggName, "_", aggLevel)
  return(assignmentResultsForDis)
}

# Computes cluster assignments for all aggregated dataseta in the directories of
# "newDataBasepathList" based on existing clustering represented by "oldDataBasePathList"
# and "clusterAssignmentsTable". For each of the new datasets (which are loaded
# from the hard-drive, need to be saved as .rds with a certain naming scheme; see
# AggregationUtility.R for the pre-aggregation routine), it is necessary to compute
# the distance the corresponding old datasets first and then assign clusters based
# on majority vote (standard is "numNeighbors = 1", so each new time series gets
# the cluster label of the closest object from the old dataset).
# The result is similar to the "analyzeAggregationImpact...()" routines, a list
# with a cluster assignment table. If "dissMatrixPath" is supplied, then the
# dissimilarity matrices are saved to the hard drive (necessary for forecasting).
assignDatasetsToClusters <- function(newDataBasePathList, oldDataBasePathList,
    clusterAssignmentsTable, disList, numNeighbors = 1, dissMatrixPath = NULL,
    expRegistry = NULL, cpuCores = 8) {
  datasetsOverviewTable <- checkAssignmentAssumptions(newDataBasePathList,
      oldDataBasePathList, clusterAssignmentsTable, disList)
  # Convert columns to avoid problems during instance generation (where we use this
  # values to access the clusterAssignmentsTable, so character is better than
  # factor, because the assignments table might have more factor levels for these
  # attributes and a factor-factor comparison might result in errors)
  datasetsOverviewTable[, Dataset := as.character(Dataset)]
  datasetsOverviewTable[, Aggregation := as.character(Aggregation)]
  setorder(datasetsOverviewTable, -Level)
  if (!is.null(dissMatrixPath)) {
    datasetsOverviewTable[, dissMatrixPath := dissMatrixPath]
    dir.create(dissMatrixPath, showWarnings = FALSE, recursive = TRUE)
  }
  if (is.null(expRegistry)) {
    expRegistry <- makeExperimentRegistry(file.dir = NA, seed = 1)
  }
  expRegistry$cluster.functions <- makeClusterFunctionsSocket(ncpus = cpuCores)
  expRegistry$source <- c("ClusteringUtility.R", "DissimilarityUtility.R")
  # Problems: different datasets (old = already clustered and new = to be assigned),
  # aggregations, levels
  addProblem(
    name = "AssignmentProblem",
    data = list(disList = disList, assignmentsTable = clusterAssignmentsTable,
                numNeighbors = numNeighbors),
    fun = function(data, job, ...) {
      return(list(
        newDataset = readRDS(file = job$pars$prob.pars$NewPath),
        oldDataset = readRDS(file = job$pars$prob.pars$OldPath),
        assignmentsTable = data$assignmentsTable[Dataset == job$pars$prob.pars$Dataset &
            Aggregation == job$pars$prob.pars$Aggregation & Level == job$pars$prob.pars$Level]
      ))
    },
    reg = expRegistry
  )
  # Algorithms: assignment for different dissimilarites and clusterings
  addAlgorithm(name = "TSClusterAssignment", fun = function(data, job, instance, ...) {
    # combine assignments for differents dissimilarties, flatten out one dimension (as
    # each result per dissimilarity itself is a named list, we simply want to concat them)
    return(unlist(lapply(instance$assignmentsTable[, as.character(unique(Dissimilarity))],
                         function(disName) {
      # Retrieve dissimilarity; name of the whole element in the dissList might be
      # different from "name" component of the element, therefore we use this code
      # to be compliant with the naming in the clustering routine; else, we could
      # also use "data$disList[[disName]]"
      relevantDisEntry <- data$disList[sapply(data$disList, function(x)
        x$name == disName)][[1]]
      # Convert assignments table -> matrix -> list of assignment vectors
      relevantAssignmentsTable <- instance$assignmentsTable[Dissimilarity == disName]
      oldAssignments <- as.matrix(relevantAssignmentsTable[, mget(grep("^A[0-9]+",
          colnames(relevantAssignmentsTable), value = TRUE))])
      oldAssignments <- lapply(seq_len(nrow(oldAssignments)), function(i)
        oldAssignments[i, ])
      names(oldAssignments) <- relevantAssignmentsTable[, Clustering]
      # Assignment routine
      assignOneDatasetToClusters(
        oldDataset = instance$oldDataset,
        newDataset = instance$newDataset,
        oldAssignmentsList = oldAssignments,
        numNeighbors = data$numNeighbors,
        datasetName = job$pars$prob.pars$Dataset,
        disEntry = relevantDisEntry,
        aggName = job$pars$prob.pars$Aggregation,
        aggLevel = job$pars$prob.pars$Level,
        dissMatrixPath = job$pars$prob.pars$dissMatrixPath
      )
    }), recursive = FALSE))
  }, reg = expRegistry)
  # Experiment: Combine problems + algorithms
  addExperiments(
    prob.designs = list(AssignmentProblem = datasetsOverviewTable),
    algo.designs = NULL,
    reg = expRegistry
  )
  submitJobs(reg = expRegistry)
  waitForJobs(reg = expRegistry)
  # Results can be retrieved as list with one element per job, but each of these
  # elements summarizes multiple cluster assignments, therefore "unlist"
  return(rbindlist(unlist(reduceResultsList(reg = expRegistry), recursive = FALSE)))
}

# Method doing the same as "assignDatasetsToClusters()". However, it does not
# parallelize over Dataset-Aggregation-Level, but Dataset-Aggregation-Level-Dissimilarity.
# If you have slow (like O(n^2)) dissimilarities, many CPU core and big datasets,
# this additional parallelization dimension might help to achieve a speed-up. It
# introduces an additional overhead, as the same dataset (with aggregation and
# level) file is read multiple times, once for each dissimilarity (while the
# original assignment routine reads it once, but applies the dissimilarities in
# an unparallelized loop). So for small datasets, the original routine might be
# faster.
assignDatasetsToClusters.dissParallelized <- function(newDataBasePathList, oldDataBasePathList,
    clusterAssignmentsTable, disList, numNeighbors = 1, dissMatrixPath = NULL,
    expRegistry = NULL, cpuCores = 8) {
  datasetsOverviewTable <- checkAssignmentAssumptions(newDataBasePathList,
      oldDataBasePathList, clusterAssignmentsTable, disList)
  # Convert columns to avoid problems during instance generation (where we use this
  # values to access the clusterAssignmentsTable, so character is better than
  # factor, because the assignments table might have more factor levels for these
  # attributes and a factor-factor comparison might result in errors)
  datasetsOverviewTable[, Dataset := as.character(Dataset)]
  datasetsOverviewTable[, Aggregation := as.character(Aggregation)]
  dissimilarityTable <- clusterAssignmentsTable[, 1,
      by = .(Dataset, Aggregation, Level, Dissimilarity)][, -"V1"]
  dissimilarityTable[, Dataset := as.character(Dataset)]
  dissimilarityTable[, Aggregation := as.character(Aggregation)]
  dissimilarityTable[, Dissimilarity := as.character(Dissimilarity)]
  datasetsOverviewTable <- merge(datasetsOverviewTable, dissimilarityTable)
  setorder(datasetsOverviewTable, -Level)
  if (!is.null(dissMatrixPath)) {
    datasetsOverviewTable[, dissMatrixPath := dissMatrixPath]
    dir.create(dissMatrixPath, showWarnings = FALSE, recursive = TRUE)
  }
  if (is.null(expRegistry)) {
    expRegistry <- makeExperimentRegistry(file.dir = NA, seed = 1)
  }
  expRegistry$cluster.functions <- makeClusterFunctionsSocket(ncpus = cpuCores)
  expRegistry$source <- c("ClusteringUtility.R", "DissimilarityUtility.R")
  # Problems: different datasets (old = already clustered and new = to be assigned),
  # aggregations, levels
  addProblem(
    name = "AssignmentProblem",
    data = list(disList = disList, assignmentsTable = clusterAssignmentsTable,
                numNeighbors = numNeighbors),
    fun = function(data, job, ...) {
      return(list(
        newDataset = readRDS(file = job$pars$prob.pars$NewPath),
        oldDataset = readRDS(file = job$pars$prob.pars$OldPath),
        assignmentsTable = data$assignmentsTable[Dataset == job$pars$prob.pars$Dataset &
            Aggregation == job$pars$prob.pars$Aggregation & Level == job$pars$prob.pars$Level &
            Dissimilarity == job$pars$prob.pars$Dissimilarity]
      ))
    },
    reg = expRegistry
  )
  # Algorithms: assignment for different dissimilarites and clusterings
  addAlgorithm(name = "TSClusterAssignment", fun = function(data, job, instance, ...) {
    # Retrieve dissimilarity; name of the whole element in the dissList might be
    # different from "name" component of the element, therefore we use this code
    # to be compliant with the naming in the clustering routine; else, we could
    # also use "data$disList[[job$pars$prob.pars$Dissimilarity]]"
    relevantDisEntry <- data$disList[sapply(data$disList, function(x)
      x$name == job$pars$prob.pars$Dissimilarity)][[1]]
    # Convert assignments table -> matrix -> list of assignment vectors
    oldAssignments <- as.matrix(instance$assignmentsTable[, mget(grep("^A[0-9]+",
        colnames(instance$assignmentsTable), value = TRUE))])
    oldAssignments <- lapply(seq_len(nrow(oldAssignments)), function(i)
      oldAssignments[i, ])
    names(oldAssignments) <- instance$assignmentsTable[, Clustering]
    # Assignment routine
    assignOneDatasetToClusters(
      oldDataset = instance$oldDataset,
      newDataset = instance$newDataset,
      oldAssignmentsList = oldAssignments,
      numNeighbors = data$numNeighbors,
      datasetName = job$pars$prob.pars$Dataset,
      disEntry = relevantDisEntry,
      aggName = job$pars$prob.pars$Aggregation,
      aggLevel = job$pars$prob.pars$Level,
      dissMatrixPath = job$pars$prob.pars$dissMatrixPath
    )
  }, reg = expRegistry)
  # Experiment: Combine problems + algorithms
  addExperiments(
    prob.designs = list(AssignmentProblem = datasetsOverviewTable),
    algo.designs = NULL,
    reg = expRegistry
  )
  submitJobs(reg = expRegistry)
  waitForJobs(reg = expRegistry)
  # Results can be retrieved as list with one element per job, but each of these
  # elements summarizes multiple cluster assignments, therefore "unlist"
  return(rbindlist(unlist(reduceResultsList(reg = expRegistry), recursive = FALSE)))
}
