# Script to merge summer time and winter time files (for 2 and 6 hour aggregates),
# taking the parts from each file which are aligned to midnight (local time)
# usage: "Rscript batchscripts/MergeSummerWinterCSVScript.R <<dataPath>"
library(data.table)

### Global params ###

dataPath <- commandArgs(trailingOnly = TRUE)[1]
outputPath <- paste0(dataPath, "processed/")

### Main ###

cat("---Started processing winter and summer files---\n")

dir.create(outputPath, showWarnings = FALSE)
fileNames <- list.files(dataPath, pattern = "\\.csv$")
aggregateNames <- unique(gsub("_winter|_summer|\\.csv", "", fileNames))
for (aggregateName in aggregateNames) {
  cat("Aggregate:", aggregateName)
  matchingFiles <- grep(aggregateName, fileNames, value = TRUE)
  if (length(matchingFiles) == 1) {
    cat(" - only one file, nothing to do.\n")
  } else if (length(matchingFiles) == 2) {# summer + winter
    winterFileName <- grep("winter", matchingFiles, value = TRUE)
    winterData <- data.table(read.csv(paste0(dataPath, winterFileName)))
    winterData[, start_hour := as.integer(substr(date_start, 12, 13))]
    winterData[, end_hour := as.integer(substr(date_end, 12, 13))]
    summerFileName <- grep("summer", matchingFiles, value = TRUE)
    summerData <- data.table(read.csv(paste0(dataPath, summerFileName)))
    summerData[, start_hour := as.integer(substr(date_start, 12, 13))]
    summerData[, end_hour := as.integer(substr(date_end, 12, 13))]
    # 2 hour and 6 hour aggregates affected by daily-saving time, 1 hour shift
    sensorData <- rbind(winterData[start_hour %% 2 == 0 & end_hour %% 2 == 0],
                        summerData[start_hour %% 2 == 0 & end_hour %% 2 == 0])
    sensorData[, c("start_hour", "end_hour") := NULL]
    setorder(sensorData, machine_name) # retain original order (within machine ordered by time)
    # Creation of file object necessary to get pure LF endings on Windows
    outputFile <- file(description = paste0(outputPath, aggregateName, ".csv"), open = "wb")
    write.table(sensorData, file = outputFile, quote = FALSE, sep = ",",
                row.names = FALSE, eol = "\n")
    close(outputFile)
    cat(" - merged.\n")
  } else {
    stop("Different naming convention assumed in script.")
  }
}

cat("---Finished processing winter and summer files---\n")
