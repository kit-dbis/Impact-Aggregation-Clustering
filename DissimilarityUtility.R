library(FastTSDistances)
library(pdc)

## Dissimilarities with pre-configured parameters
# each entry consists of (symbolic) name; method name; parameters; functions and
# packages (if needed for export to run parallel computations); flag if only for
# pairs of time series or whole list
# - Complexity-invariant Distance (combined with DTW, L2)
# - Compression-based dissimilarity (CDM)
# - (Temporal) Correlation-corrected Distance (combined with DTW, L2)
# - Dynamic time warping (with/without window)
# - Edit Distance with Real Penalty
# - L norms (Manhattan, Euclidean, Chebyscheff)
# - Permutation Distribution Distance
# - Shaped-Based Distance (uses cross-correlation)
DEFAULT_DISSIMILARITIES <- list(
  L1 = list(
    name = "L1", method = "l1Dist_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  L2 = list(
    name = "L2", method = "l2Dist_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  L2.CID = list(
    name = "L2.CID", method = "l2Dist_fast", params = list(cid = TRUE),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  L2.CORT = list(
    name = "L2.CORT", method = "l2Dist_fast", params = list(cortK = 2),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  Lmax = list(
    name = "Lmax", method = "lmaxDist_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  DTW = list(
    name = "DTW", method = "DTWDist_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  DTW.CID = list(
    name = "DTW.CID", method = "DTWDist_fast", params = list(cid = TRUE),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  DTW.CORT = list(
    name = "DTW.CORT", method = "DTWDist_fast", params = list(cortK = 2),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  DTW.Band = list(
    name = "DTW.Band", method = "DTWDistSakoeChiba_fast", params = list(windowSize = .1),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  ERP = list(
    name = "ERP", method = "ERPDist_fast", params = list(g = 0),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  SBD = list(
    name = "SBD", method = "shapeBasedDistance", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  CDM = list(
    name = "CDM", method = "compDistTSList",  params = list(symbolCount = 8),
    exportFunctions = c(), exportPackages = c(), takesList = TRUE
  ),
  PDD = list(
    name = "PDD", method = "pdcDistTSList",  params = list(),
    exportFunctions = c(), exportPackages = c("pdc"), takesList = TRUE
  )
)

## Multivariate issimilarities with pre-configured parameters
# same dissimilarities and same structure as above
DEFAULT_MULTIVARIATE_DISSIMILARITIES <- list(
  L1 = list(
    name = "L1", method = "l1DistMult_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  L2 = list(
    name = "L2", method = "l2DistMult_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  L2.CID = list(
    name = "L2.CID", method = "l2DistMult_fast", params = list(cid = TRUE),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  L2.CORT = list(
    name = "L2.CORT", method = "l2DistMult_fast", params = list(cortK = 2),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  Lmax = list(
    name = "Lmax", method = "lmaxDistMult_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  DTW = list(
    name = "DTW", method = "DTWDistMult_fast", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  DTW.CID = list(
    name = "DTW.CID", method = "DTWDistMult_fast", params = list(cid = TRUE),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  DTW.CORT = list(
    name = "DTW.CORT", method = "DTWDistMult_fast", params = list(cortK = 2),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  DTW.Band = list(
    name = "DTW.Band", method = "DTWDistSakoeChibaMult_fast", params = list(windowSize = .1),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  ERP = list(
    name = "ERP", method = "ERPDistMult_fast", params = list(g = 0),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  SBD = list( #same function for uni- and multivariate
    name = "SBD", method = "shapeBasedDistance", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  CDM = list( #same function for uni- and multivariate
    name = "CDM", method = "compDistTSList",  params = list(symbolCount = 8),
    exportFunctions = c(), exportPackages = c(), takesList = TRUE
  ),
  PDD = list(
    name = "PDD", method = "pdcDistTSListMult",  params = list(),
    exportFunctions = c(), exportPackages = c("pdc"), takesList = TRUE
  )
)

## Dissimilarities with pre-configured parameters which can handle univariate and
# multivariate time series
# same dissimilarities and same structure as above
DEFAULT_DISPATCH_DISSIMILARITIES <- list(
  L1 = list(
    name = "L1", method = "l1DistDispatcher", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  L2 = list(
    name = "L2", method = "l2DistDispatcher", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  L2.CID = list(
    name = "L2.CID", method = "l2DistDispatcher", params = list(cid = TRUE),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  L2.CORT = list(
    name = "L2.CORT", method = "l2DistDispatcher", params = list(cortK = 2),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  Lmax = list(
    name = "Lmax", method = "lmaxDistDispatcher", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  DTW = list(
    name = "DTW", method = "DTWDistDispatcher", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  DTW.CID = list(
    name = "DTW.CID", method = "DTWDistDispatcher", params = list(cid = TRUE),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  DTW.CORT = list(
    name = "DTW.CORT", method = "DTWDistDispatcher", params = list(cortK = 2),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  DTW.Band = list(
    name = "DTW.Band", method = "DTWDistSakoeChibaDispatcher", params = list(windowSize = .1),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  ERP = list(
    name = "ERP", method = "ERPDistDispatcher", params = list(g = 0),
    exportFunctions = c(), exportPackages = c(), takesList = FALSE
  ),
  SBD = list( #same function for uni- and multivariate
    name = "SBD", method = "shapeBasedDistance", params = list(), exportFunctions = c(),
    exportPackages = c(), takesList = FALSE
  ),
  CDM = list( #same function for uni- and multivariate
    name = "CDM", method = "compDistTSList",  params = list(symbolCount = 8),
    exportFunctions = c(), exportPackages = c(), takesList = TRUE
  ),
  PDD = list(
    name = "PDD", method = "pdcDistTSListDispatcher",  params = list(),
    exportFunctions = c(), exportPackages = c("pdc"), takesList = TRUE
  )
)

### Functions ###

# Dispatching functions which support univariate as well as multivariate time series

DTWDistDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    DTWDistMult_fast(x, y, ...)
  } else {
    DTWDist_fast(x, y, ...)
  }
}

DTWDistSakoeChibaDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    DTWDistSakoeChibaMult_fast(x, y, ...)
  } else {
    DTWDistSakoeChiba_fast(x, y, ...)
  }
}

EDRDistDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    EDRDistMult_fast(x, y, ...)
  } else {
    EDRDist_fast(x, y, ...)
  }
}

ERPDistDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    ERPDistMult_fast(x, y, ...)
  } else {
    ERPDist_fast(x, y, ...)
  }
}

l1DistDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    l1DistMult_fast(x, y, ...)
  } else {
    l1Dist_fast(x, y, ...)
  }
}

l2DistDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    l2DistMult_fast(x, y, ...)
  } else {
    l2Dist_fast(x, y, ...)
  }
}

lmaxDistDispatcher <- function(x, y, ...) {
  if (is.matrix(x)) {
    lmaxDistMult_fast(x, y, ...)
  } else {
    lmaxDist_fast(x, y, ...)
  }
}

pdcDistTSListDispatcher <- function(tsList, ...) {
  # method can be called with a list of time series or a list of two lists of time series
  if (is.list(tsList[[1]])) {
    if (is.matrix(tsList[[1]][[1]])) {
      pdcDistTSListMult(tsList, ...)
    } else {
      pdcDistTSList(tsList, ...)
    }
  } else if (is.matrix(tsList[[1]])) {
    pdcDistTSListMult(tsList, ...)
  } else {
    pdcDistTSList(tsList, ...)
  }
}
