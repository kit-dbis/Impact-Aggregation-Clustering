# Call with "Rscript batchscripts/RandomDataGenerationScript.R <<arguments>>"
# from project directory
# Supply the following arguments:
# --dataPath="path/to/dir/" : directory where the generated data should be saved
# --days=60 : number of days
# --machines=10 : number of machines
# --samplingInterval=3 : interval between measurements in the "original" data in secs

library(data.table)
library(FastTSDistances)
library(optparse)

source("PreprocessingUtility.R")

## Reproducability
set.seed(25)

## Input params
params <- parse_args(object = OptionParser(option_list = list(
    make_option("--dataPath", action = "store", type = "character", default = "./"),
    make_option("--days", action = "store", type = "integer", default = 60),
    make_option("--machines", action = "store", type = "integer", default = 10),
    make_option("--samplingInterval", action = "store", type = "numeric", default = 3)
  ), add_help_option = FALSE))

## Create un-aggregated dataset table (could also be fed into aggregation
# analysis directly, triggering automatic computation of aggregates; but, as in
# our other experiments, we will not use this original dataset, but only aggregates
# which are created in the 2nd part of this script)
cat("Creating un-aggregated data.\n")
numOfMeasurements <- params$days * 60*60*24 / params$samplingInterval
# machines are named AVT_01 etc.; add leading zeros depending on number of machines;
# use "each" parameter so all timestamps will be added later for all machines automatically
randomData <- data.table(machine_name = as.factor(rep(sprintf(paste0("AVT_%0",
    nchar(params$machines), "d"), seq_len(params$machines)), each = numOfMeasurements)))
invisible(gc()) # without "gc()", memory is not released fast enough in our sequence of computations
randomData[, sensor_name := factor("DummySensor")] # same sensor for all
# Add timestamps; these are replicated automatically for all machines
randomData[, c("ts", "ts_numeric", "day", "secOfDay") := createTimestamps.aggData(
  from = as.POSIXct("2017-01-01", tz = "GMT"), samplingInterval = params$samplingInterval,
  numberOfTS = numOfMeasurements
)]
invisible(gc())
stopifnot(all(randomData[, uniqueN(ts_numeric), by = machine_name]$V1 == numOfMeasurements))
randomData[, mean := rnorm(n = .N)] # not really a mean, but for compatibility
invisible(gc())
setorder(randomData, ts_numeric, machine_name)
invisible(gc())
generatedDataPath <- paste0(params$dataPath, "generatedData/")
dir.create(generatedDataPath, showWarnings = F, recursive = T)
cat("Saving un-aggregated data.\n")
saveRDS(randomData, paste0(generatedDataPath, "randomData.rds"))

## Create aggregate files (similar to our real-world data, aggregates are computed
# on even finer measurements [here: randomData], but only aggregates used in experiments)
# intervals are same as for real-world data: 30 s, 1/5/10/15/30 min, 1/2/6 h
for (aggInterval in c(30, 60, 300, 600, 900, 1800, 3600, 7200, 21600)) {
  cat("Aggregating to ", aggInterval, " second intervals.\n", sep = "")
  aggData <- randomData[, .(ts = ts[1], ts_numeric = ts_numeric[1], secOfDay = secOfDay[1],
      median = PMedAA_fast(mean, 1), mean = PAA_fast(mean, 1), min = PMinAA_fast(mean, 1),
      max = PMaxAA_fast(mean, 1), stddev = PSDAA_fast(mean, 1, sample = T),
      skewness = PSkewAA_fast(mean, 1), kurtosis = PKurtAA_fast(mean, 1, excess = T),
      count = .N, sample = sample(mean, size = 1)),
    by = .(machine_name, sensor_name, day, V1 = floor(secOfDay / aggInterval))]
  aggData[, random := rnorm(n = .N)]
  aggData[, V1 := NULL]
  gc()
  setorder(randomData, ts_numeric, machine_name)
  gc()
  saveRDS(aggData, paste0(params$dataPath, "RandomDataAgg_", aggInterval, "_sec.rds"))
  gc()
}
