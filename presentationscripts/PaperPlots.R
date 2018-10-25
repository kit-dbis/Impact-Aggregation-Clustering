# Not intended for general use, contains hard-coded paths and is tied to certain datatsets.
source("BatchExperimentEvaluationUtility.R")
source("PreprocessingUtility.R")
source("presentationscripts/PaperPlotUtility.R")

library(data.table)
library(ggplot2)

### Global params ###

DATA_BASE_PATH <- "<<path/on/your/PC>>/clustagg-data/data/"
PLOT_OUT_PATH <- "<<path/on/your/PC>>"
# Store merged experimental results
FLENGTH_DATA_OUT_PATH <- paste0(DATA_BASE_PATH, "MergedFixedLength/")
VLENGTH_DATA_OUT_PATH <- paste0(DATA_BASE_PATH, "MergedVariableLength/")
MERGED_DATA_OUT_PATH <- paste0(DATA_BASE_PATH, "Merged/")

### Prepare merged datasets used for evaluation and load ###

RANDOM_DATA_DIR <- paste0(DATA_BASE_PATH, "Gaussian600TS/")
FLENGTH_DATA_DIRS <- list.files(DATA_BASE_PATH, pattern = "^2017_Jan_Aug_", full.names = TRUE)
FLENGTH_DATA_DIRS <- addSlashesToDirs(FLENGTH_DATA_DIRS[!grepl(pattern = "VariableLength", FLENGTH_DATA_DIRS)])
stopifnot(length(FLENGTH_DATA_DIRS) == 6) # six physical quantities used in experiments
VLENGTH_DATA_DIRS <- addSlashesToDirs(list.files(DATA_BASE_PATH, pattern = "VariableLength", full.names = TRUE))
stopifnot(length(VLENGTH_DATA_DIRS) == 6) # six physical quantities used in experiments

preparePaperPlotDataset(fLengthDirs = FLENGTH_DATA_DIRS, vLengthDirs = VLENGTH_DATA_DIRS,
    fLengthOutDir = FLENGTH_DATA_OUT_PATH, vLengthOutDir = VLENGTH_DATA_OUT_PATH,
    mergedOutDir = MERGED_DATA_OUT_PATH_PATH) # merge, rename some experimental settings
for (resultsFileName in list.files(paste0(MERGED_DATA_OUT_PATH_PATH, "results/"),
                                   full.names = TRUE, pattern = "\\.rds")) {
  assign(gsub(pattern = ".rds", replacement = "", x = basename(resultsFileName), fixed = TRUE),
         value = readRDS(resultsFileName))
}

###############################################################################

##### Experimental Design #####

### Table I: Overview of datasets ###

summarizeTimeSeriesCount(FLENGTH_DATA_DIRS, print = FALSE)
summarizeTimeSeriesCount(VLENGTH_DATA_DIRS, print = FALSE)

##### Evaluating the Effects of Aggregation #####

### Figure 1: Validity evaluation approaches ###
# (created with IPE)
brewer.pal(6, "PRGn")[1] # purple
brewer.pal(6, "PRGn")[5] # green
brewer.pal(6, "Set3")[6] # orange

##### Results #####

### Figure 2(a): Impact on internal and external (base) by aggregation ###

plotTable <- rbind(internalCVITable[Index == "Sil" & Reference == "B" &
      Level != 2880, .(Value = median(Value, na.rm = TRUE)), by = .(Aggregation, Level, Index)],
    externalCVITable[Index == "i.v-D" & Reference == "B" &
      Level != 2880, .(Value = median(Value, na.rm = TRUE)), by = .(Aggregation, Level, Index)])
plotTable[Index == "Sil", Index := "Internal (Silhouette)"]
plotTable[Index == "i.v-D", Index := "External (Inv. van Dongen)"]
impactLinePlotForPaper(plotTable, yAxisName = "Base level validity") +
  theme(axis.title.y = element_text(hjust = 0.7), legend.margin = margin(0,100,0,0)) +
  labs(col = "Aggre-\ngation", shape = "Aggre-\ngation", linetype = "Aggre-\ngation")
saveImpactLinePlot(paste0(PLOT_OUT_PATH, "ImpactBaseByAggregation.pdf"))

### Figure 2(b): Impact on external (previous) by dissimilarity ###

impactLinePlotForPaper(externalCVITable[Index == "i.v-D" & Reference == "P" &
    Level != 2880 & !Dissimilarity %in% c("DTW.Band", "DTW.CID", "DTW.CORT", "L2.CID", "L2.CORT"),
    .(Value = median(Value, na.rm = TRUE)), by = .(Dissimilarity, Level, vLength)],
    yAxisName = "Inv. van Dongen") +
  theme(axis.title.y = element_text(hjust = 0.7))
saveImpactLinePlot(paste0(PLOT_OUT_PATH, "ImpactExternalPreviousByDissimilarity.pdf"))

### Not a figure any more: Correlation between indices ###

# Internal (remove NAs and infinity values)
corrMatrix <- summarizeExperimentSettingCorrelation(internalCVITable[Index != "D-B"],
    category = "Index", plot = FALSE, infStrategy = "drop", correlationMethod = "spearman")
corrPlotForPaper(corrMatrix)
saveCorrMatrixPlot(paste0(PLOT_OUT_PATH, "InternalIndicesCorrelation.pdf"))
# External (remove NAs, exclude highest level)
corrMatrix <- summarizeExperimentSettingCorrelation(externalCVITable[Reference %in% c("B", "P") & Level != 2880],
  category = "Index", plot = FALSE, infStrategy = "drop", correlationMethod = "spearman")
corrPlotForPaper(corrMatrix)
saveCorrMatrixPlot(paste0(PLOT_OUT_PATH, "ExternalIndicesCorrelation.pdf"))

### Figure 3(a): Impact on internal (current) as boxplot ###

plotTable <- internalCVITable[Index %in% c("Sil", "i.Con") &
    Reference == "C" & !is.na(Value), .(Level, Value, Index)]
plotTable[Index == "Sil", Index := "Silhouette"]
plotTable[Index == "i.Con", Index := "Inverted Connectivity"]
impactBoxPlotForPaper(plotTable, yAxisName = "Internal validity") +
  theme(axis.title.y = element_text(hjust = 1))
saveImpactBoxPlot(paste0(PLOT_OUT_PATH, "ImpactBoxInternalCurrent.pdf"))

### Figure 3(b): Impact on internal (base) as boxplot ###

plotTable <- internalCVITable[Index %in% c("Sil", "i.Con") &
    Reference == "B" & !is.na(Value), .(Level, Value, Index)]
plotTable[Index == "Sil", Index := "Silhouette"]
plotTable[Index == "i.Con", Index := "Inverted Connectivity"]
impactBoxPlotForPaper(plotTable, yAxisName = "Internal validity") +
  theme(axis.title.y = element_text(hjust = 1))
saveImpactBoxPlot(paste0(PLOT_OUT_PATH, "ImpactBoxInternalBase.pdf"))

### Figure 4: Impact on internal (current) for random data ###

randomDataInternalCVITable <- readRDS(paste0(RANDOM_DATA_DIR, "results/internalCVITable.rds"))
plotTable <- randomDataInternalCVITable[Index == "Silhouette" & Reference == "C" &
    Aggregation %in% c("mean", "stddev", "random", "sample"),
    .(Value = median(Value, na.rm = TRUE)), by = .(Aggregation, Level, Dataset)]
plotTable <- rbind(plotTable, internalCVITable[Index == "Sil" & Reference == "C" &
    Aggregation %in% c("mean", "stddev", "random", "sample"),
    .(Value = median(Value, na.rm = TRUE)), by = .(Aggregation, Level)], fill = TRUE)
plotTable[Dataset == "GaussianRandomData", Dataset := "Gaussian random data"]
plotTable[is.na(Dataset), Dataset := "Energy data"]
impactLinePlotForPaper(plotTable, yAxisName = "Silhouette Coefficient") +
  theme(axis.title.y = element_text(hjust = 0.7), legend.margin = margin(0,80,0,0))
saveImpactLinePlot(paste0(PLOT_OUT_PATH, "ImpactInternalCurrentRandomData.pdf"))
