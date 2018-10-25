library(data.table)
library(extrafont)
# extrafont::font_import() # needs only to be done once after package installation
library(ggplot2)
library(RColorBrewer)

source("BatchExperimentEvaluationUtility.R")

### Global params ###

Sys.setenv(R_GSCMD = "<<path/to/Ghostscript>>/bin/gswin64c.exe") # for font embedding

### Functions ###

# Saves a plot with standard size and embeds fonts
savePaperPlot <- function(plotPath, plotWidth = 8, plotHeight = 4.5,
                           plotDevice = "pdf") {
  ggsave(file = paste0(plotPath), width = plotWidth, height = plotHeight,
         device = plotDevice)
  if (plotDevice == "pdf") {
    embed_fonts(file = plotPath)
  }
  return(TRUE)
}

# Saves a plot in square format and embeds fonts
saveCorrMatrixPlot <- function(plotPath, plotSize = 6, plotDevice = "pdf") {
  savePaperPlot(plotPath, plotWidth = plotSize, plotHeight = plotSize,
                 plotDevice = plotDevice)
}

# Saves an impact line plot in appropriate format and embeds fonts
saveImpactLinePlot <- function(plotPath, plotWidth = 12, plotHeight = 6.5,
    plotDevice = "pdf") {
  savePaperPlot(plotPath, plotWidth = plotWidth, plotHeight = plotHeight,
                 plotDevice = plotDevice)
}

# Saves an impact box plot in appropriate format and embeds fonts
saveImpactBoxPlot <- function(plotPath, plotWidth = 12, plotHeight = 5.5,
    plotDevice = "pdf") {
  savePaperPlot(plotPath, plotWidth = plotWidth, plotHeight = plotHeight,
                plotDevice = plotDevice)
}

# Adds a default plot styling for plots used in the paper
paperPlotTheme <- function(evaluationTable = NULL, fontSize = 28) {
  result <- list(
    scale_colour_brewer(palette = "Dark2"),
    scale_fill_brewer(palette = "Dark2"),
    theme_bw(), # black-and-white theme
    theme(
      text = element_text(size = fontSize, family = "Linux Biolinum G"),
      axis.text = element_text(size = fontSize),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = fontSize),
      legend.key.size = unit(fontSize, "pt"))) # legend symbols
  if (!is.null(evaluationTable) && evaluationTable[!is.na(Value), all(Value >= 0 & Value <= 1)]) {
    result <- c(result, list(scale_y_continuous(limits = c(0,1),
        breaks = function(x) seq(from = 0, to = 1, by = 0.2),
        minor_breaks = function(x) seq(from = 0, to = 1, by = 0.1))))
  }
  return(result)
}

# based on code from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Creates a correlation plot
corrPlotForPaper <- function(corrMatrix, corrNumberSize = 10, axisFontSize = 50,
    style = paperCorrPlotTheme(axisFontSize = axisFontSize)) {
  corrMatrix <- round(corrMatrix, 2)
  corrTable <- data.table(melt(corrMatrix, na.rm = TRUE))
  corrTable[, Var1 := as.factor(Var1)]
  corrTable[, Var2 := as.factor(Var2)]
  paletteColors <- brewer.pal(3, "Dark2")
  ggplot(data = corrTable, mapping = aes(x = Var1, y = Var2, fill = value, label = value)) +
    scale_fill_gradient2(low = paletteColors[2], high = paletteColors[1], mid = "white", 
                         midpoint = 0, limit = c(-1,1)) +
    coord_fixed() + # squares
    geom_tile(color = "white") + # background grid between squares
    geom_text(size = corrNumberSize, fontface = "bold") +
    style
}

# Theme for correlation matrix plots
paperCorrPlotTheme <- function(axisFontSize) {
  return(list(theme_minimal(), # no gray background
      theme(
        text = element_text(size = axisFontSize, family = "Linux Biolinum G"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
  ))
}

# Plots all evaluation index values from an "evaluationTable" (produced by
# "createClusteringEvaluationTable()") against the aggregation level;
# adaptation of ClusteringEvaluationUtility::plotClusteringEvaluationOverview()
# which does some extra formatting;
# requires at least a "Value" and a "Level" column; if a third column exists, its
# values will be plotted as group (by color)
impactBoxPlotForPaper <- function(evaluationTable, yAxisName = "Value", lineSize = 1,
    fontSize = 40, style = paperImpactBoxPlotTheme(fontSize = fontSize, evaluationTable = evaluationTable)) {
  evaluationTable <- convertLevelToSamplingIntervals(evaluationTable)
  groupCol <- setdiff(colnames(evaluationTable), c("Value", "Level"))
  if (!("Level" %in% colnames(evaluationTable)) || !("Value" %in% colnames(evaluationTable)) ||
      length(groupCol) > 1) {
    stop(paste0("Routine requires 2 or 3 columns, \"Level\", \"Value\" and optionally ",
                "a third one used for grouping"))
  }
  if (length(groupCol) == 1) {
    result <- ggplot(data = evaluationTable, mapping = aes(x = Level, y = Value,
        color = get(groupCol))) + labs(color = groupCol)
    result <- result + geom_boxplot(size = lineSize, outlier.size = lineSize) + style +
      scale_color_manual(values = c(brewer.pal(6, "PRGn")[1], brewer.pal(6, "PRGn")[5]))
  } else {# No group column
    result <- ggplot(data = evaluationTable, mapping = aes(x = Level, y = Value))
    result <- result + geom_boxplot(color = brewer.pal(3, "Dark2")[2], size = lineSize,
                                    outlier.size = lineSize) + style
  }
  result <- result + labs(x = "Aggregation level", y = yAxisName)
  return(result)
}

# Theme for impact box plots
paperImpactBoxPlotTheme <- function(evaluationTable, fontSize) {
  return(c(paperPlotTheme(evaluationTable, fontSize = fontSize), list(theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))))
}

# Replace numeric Levels with a factor whose labels represent sampling intervals
convertLevelToSamplingIntervals <- function(clusterSummaryTable) {
  clusterSummaryTable <- copy(clusterSummaryTable)
  samplingIntervals <- c("2880" = "30 s", "1440" = "1 min", "288" = "5 min",
      "144" = "10 min", "96" = "15 min", "48" = "30 min", "24" = "1 h",
      "12" = "2 h", "4" = "6 h")
  aggregationLevels <- clusterSummaryTable[, sort(unique(Level), decreasing = TRUE)]
  clusterSummaryTable[, Level := factor(Level, levels = aggregationLevels,
      labels = samplingIntervals[as.character(aggregationLevels)], ordered = TRUE)]
  return(clusterSummaryTable)
}

# Create a plot which shows the value of one evaluation measure over the aggregation
# levels, grouped by one experimental setting; requires an "evaluationTable" with
# the columns "Value", "Level" and a third one which represent the group;
# optionally, a forth column might be provided which allows to add a facet;
# adaptation of ClusteringEvaluationUtility::plotOneClusteringEvaluation()
# which is simplified and does some extra formatting
impactLinePlotForPaper <- function(evaluationTable, yAxisName = "Value", fontSize = 40,
    style = paperImpactLinePlotTheme(evaluationTable = evaluationTable, fontSize = fontSize)) {
  groupCol <- setdiff(colnames(evaluationTable), c("Value", "Level"))
  if (length(groupCol) > 2) {
    stop(paste0("Group column is ambiguous. Make sure your table only contains ",
                "three or four columns, two being \"Value\" and \"Level\"."))
  }
  if (evaluationTable[, .N, by = setdiff(colnames(evaluationTable), "Value")][, any(N > 1)]) {
    stop("There are several values for one point to plot.")
  }
  # Replace level numbers with sorted factor denoting sampling intervals
  evaluationTable <- convertLevelToSamplingIntervals(evaluationTable)
  # Lines
  result <- ggplot() +
     geom_line(data = evaluationTable, aes(x = Level, y = Value, color = get(groupCol[1]),
         linetype = get(groupCol[1]), group = get(groupCol[1])), size = 1) +
     geom_point(data = evaluationTable, aes(x = Level, y = Value, color = get(groupCol[1]),
         shape = get(groupCol[1])), size = 3) +
     # Allow for more than 6 shapes
     scale_shape_manual(values = seq_len(evaluationTable[, uniqueN(get(groupCol[1]))])) +
     labs(linetype = groupCol[1], shape = groupCol[1]) # map all to same legend
  # Axis stuff
  if (length(groupCol) == 2) {
    if (groupCol[2] == "vLength") {
      labelFunc <- as_labeller(c("FALSE" = "Fixed length", "TRUE" = "Variable length"))
    } else {
      labelFunc <- "label_value" # ggplot default
    }
    result <- result + facet_grid(as.formula(paste0(".~", groupCol[2])), labeller = labelFunc)
  }
  # Styling
  result <- result +
    theme(legend.position = "right", legend.spacing = unit(0, "pt")) +
    labs(x = "Aggregation level", y = yAxisName, color = groupCol[1]) +
    style
  return(result)
}

# Theme for impact line plots
paperImpactLinePlotTheme <- function(evaluationTable, fontSize) {
  return(c(paperPlotTheme(evaluationTable, fontSize = fontSize), list(theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))))
}

# Merges results for experiments with time series of fixed length ("flengthDirs")
# and variable length ("vLengthDirs"); write the output tables and then also
# merges over these categories, also renaming experimental settings for plots.
preparePaperPlotDataset <- function(fLengthDirs, vLengthDirs, fLengthOutDir,
    vLengthOutDir, mergedOutDir) {
  if (length(list.files(paste0(fLengthOutDir, "results/"), pattern = "\\.rds$")) == 0) {
    cat("### Merging results with time series of fixed length. ###\n")
    saveAllMergedEvaluationTables(dataPaths = fLengthDirs, outputPath = fLengthOutDir)
    clusterCountTable <- prepareClusterCountTableForPlots(mergedDataPath = fLengthOutDir,
        inputDataPaths = fLengthDirs) # add more columns to cluster count table
    saveRDS(clusterCountTable, file = paste0(fLengthOutDir, "results/clusterCountTable.rds"))
  }
  if (length(list.files(paste0(vLengthOutDir, "results/"), pattern = "\\.rds$")) == 0) {
    cat("### Merging results with time series of variable length. ###\n")
    saveAllMergedEvaluationTables(dataPaths = vLengthDirs, outputPath = vLengthOutDir)
    clusterCountTable <- prepareClusterCountTableForPlots(mergedDataPath = vLengthOutDir,
        inputDataPaths = vLengthDirs)
    saveRDS(clusterCountTable, file = paste0(vLengthOutDir, "results/clusterCountTable.rds"))
  }
  if (length(list.files(paste0(mergedOutDir, "results/"), pattern = "\\.rds$")) == 0) {
    cat("### Merging results with time series of fixed and variable length. ###\n")
    saveAllMergedEvaluationTables(dataPaths = c(fLengthOutDir, vLengthOutDir),
                                  outputPath = mergedOutDir)
    cat("### Rename table content for plots. ###\n")
    for (resultsFileName in list.files(paste0(mergedOutDir, "results/"),
                                       full.names = TRUE, pattern = "\\.rds")) {
      cat(gsub(pattern = ".rds", replacement = "", x = basename(resultsFileName), fixed = TRUE), "\n")
      clusterSummaryTable <- readRDS(resultsFileName)
      clusterSummaryTable <- addVLengthColumn(clusterSummaryTable)
      clusterSummaryTable <- renameForPlots(clusterSummaryTable)
      saveRDS(clusterSummaryTable, file = resultsFileName)
    }
  }
}

# Replace names in all relevant columns to make them more plot-friendly.
# Inlines the different renaming routines to reduce copy, format and sort actions.
renameForPlots <- function(clusterSummaryTable) {
  clusterSummaryTable <- copy(clusterSummaryTable)
  # Dataset
  clusterSummaryTable[Dataset %like% "ActivePower", Dataset := "Active power"]
  clusterSummaryTable[Dataset %like% "Amperage", Dataset := "Amperage"]
  clusterSummaryTable[Dataset %like% "Frequency", Dataset := "Frequency"]
  clusterSummaryTable[Dataset %like% "PowerFactor", Dataset := "Power factor"]
  clusterSummaryTable[Dataset %like% "PositiveEnergy", Dataset := "Positive energy"]
  clusterSummaryTable[Dataset %like% "Voltage", Dataset := "Voltage"]
  # Dissimilarity
  clusterSummaryTable[Dissimilarity == "CDM", Dissimilarity := "CDM+"]
  # Cluster validity
  if ("Index" %in% colnames(clusterSummaryTable)) {
    # External CVI
    clusterSummaryTable[Index == "inv.vanDongen", Index := "i.v-D"]
    clusterSummaryTable[Index == "Fowlkes-Mallows", Index := "F-M"]
    # Internal CVI
    clusterSummaryTable[Index == "Silhouette", Index := "Sil"]
    clusterSummaryTable[Index == "inv.Connectivity", Index := "i.Con"]
    clusterSummaryTable[Index == "Gen.Dunn", Index := "Dunn"]
    clusterSummaryTable[Index == "Gen.Davies.Bouldin", Index := "D-B"]
    clusterSummaryTable[Index == "inv.Gen.Davies.Bouldin", Index := "i.D-B"]
    clusterSummaryTable[ , Index := as.factor(as.character(Index))] # remove levels
  }
  # Forecasting
  if ("Reference" %in% colnames(clusterSummaryTable)) {
    clusterSummaryTable[Reference == "naive", Reference := "Naive"]
    clusterSummaryTable[Reference == "clus.weighted", Reference := "Clustering, weighted"]
    clusterSummaryTable[Reference == "clus.mean", Reference := "Clustering, unweighted"]
    clusterSummaryTable[Reference == "all.weighted", Reference := "All, weighted"]
    clusterSummaryTable[Reference == "all.mean", Reference := "All, unweighted"]
    clusterSummaryTable[Reference == "10NN.weighted", Reference := "10NN, weighted"]
    clusterSummaryTable[Reference == "10NN.mean", Reference := "10NN, unweighted"]
    clusterSummaryTable[, Reference := as.factor(as.character(Reference))] # remove levels
  }
  clusterSummaryTable <- formatAndSortClusterSummaryTable(clusterSummaryTable)
  return(clusterSummaryTable)
}

# Takes a table summarizing clustering results containing datasets of fixed and
# variable length and denotes this property by a separate column, renaming the
# datasets (so time series from the same physical quantity, fixed and variable,
# can be distinguished in plots).
addVLengthColumn <- function(clusterSummaryTable) {
  clusterSummaryTable <- copy(clusterSummaryTable)
  clusterSummaryTable[, vLength := grepl(pattern = "VariableLength", x = Dataset)]
  clusterSummaryTable[, Dataset := gsub(pattern = "VariableLength", replacement = "", x = Dataset)]
  clusterSummaryTable <- formatAndSortClusterSummaryTable(clusterSummaryTable)
  return(clusterSummaryTable)
}
