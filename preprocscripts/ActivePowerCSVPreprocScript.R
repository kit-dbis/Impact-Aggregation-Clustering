# According to a domain expert, active power data should be positive, at least if
# voltage and amperage are positive. We have some machines whose sensor values are
# completely negative, but apart from that show similar behaviour. As a consequence,
# this script inverts these sensor data, i.e., changes the sign. Machines that have
# positive as well as negative values are not affected.
library(data.table)
source("PreprocessingUtility.R")

dataPath <- "../../data/2017_Jan_Aug_ActivePower1456Minmax/"

# Standard CSV -> RDS conversion
preprocessAggDataDir(dataPath)

# Find machines with inverted value
inputRDSFiles <- list.files(path = dataPath, pattern = "\\.rds$", full.names = TRUE)
inputRDSFiles <- inputRDSFiles[!endsWith(inputRDSFiles, "_valid.rds")]
datasetFile <- inputRDSFiles[which.max(file.size(inputRDSFiles))]
dataset <- readRDS(datasetFile)
invertedMachines <- dataset[, all(mean <= 0), by = machine_name][V1 == TRUE, as.character(machine_name)]

for (inputFilePath in inputRDSFiles) {
  cat("File:", inputFilePath, "\n")
  outputFilePath <- gsub("\\.rds$", "_valid.rds", inputFilePath)
  if (file.exists(outputFilePath)) {
    cat("Already pre-processed.\n")
  } else {
    dataset <- readRDS(inputFilePath)
    invertedMachinesTable <- dataset[machine_name %in% invertedMachines]
    invertedMachinesTable[, mean := -mean]
    invertedMachinesTable[, median := -median]
    invertedMachinesTable[, min := -min]
    invertedMachinesTable[, max := -max]
    setnames(invertedMachinesTable, c("min", "max"), c("max", "min"))
    invertedMachinesTable[, skewness := -skewness]
    # standard deviation, kurtosis, count not affected
    dataset <- rbind(dataset[!(machine_name %in% invertedMachines)], invertedMachinesTable)
    setorder(dataset, ts_numeric, machine_name)
    saveRDS(dataset, outputFilePath)
  }
}
