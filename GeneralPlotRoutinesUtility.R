library(data.table)
library(doParallel)
library(foreach)
library(ggplot2)
library(tcltk)

### Functions ###

# Calls geom_point for tables with one row and geom_line else
geom_dispatch <- function(...) {
  dots <- list(...)
  if (nrow(dots$data) == 1) {
    geom_point(...)
  } else {
    geom_line(...)
  }
}

# Plots all sensors from a data.table (containing machine name, sensor name,
# timestamp, value) to a "directory"; depending on "combineSameSensors", the
# same sensor in different machines will be plotted together or separately
plotSensorsOfMachines <- function(sensorDataTable, combineSameSensors = TRUE,
    directory = "./", plotWidth = 960, plotHeight = 540) {
  dir.create(path = directory, recursive = TRUE, showWarnings = FALSE)
  sensorNames <- sensorDataTable[, unique(sensor_name)]
  sensorCount <- length(sensorNames)
  for (i in 1:sensorCount) {
    if (!exists("progBar")) {
      progBar <- tkProgressBar(max = sensorCount, title = "Plot progress")
    }
    setTkProgressBar(progBar, value = i, label = paste0("Creating plot(s) for sensor ",
                                                        i, "/", sensorCount))
    oneSensorData <- sensorDataTable[sensor_name == sensorNames[i]]
    if (combineSameSensors) {
      png(filename = paste0(directory, sensorNames[i], ".png"),
          width = plotWidth, height = plotHeight)
      print(ggplot(oneSensorData) +
              geom_line(aes(x = ts, y = value, color = machine_name)) +
              ggtitle(sensorNames[i]))
      dev.off()
    } else {
      for (machineName in oneSensorData[, unique(machine_name)]) {
        combinationName <- paste0(machineName, "_", sensorNames[i])
        png(filename = paste0(directory, combinationName, ".png"),
            width = plotWidth, height = plotHeight)
        print(ggplot(oneSensorData[machine_name == machineName]) +
                geom_line(aes(x = ts, y = value)) +
                ggtitle(combinationName))
        dev.off()
      }
    }
    gc()
  }
  close(progBar)
  return(TRUE)
}

# Plots all time series from the "timeSeriesList" to the "directory";
# if the time series items are named, the names are assumed to be POSIXct
# dates - else the items are just numbered 
plotTimeSeriesList <- function(timeSeriesList, directory = "./", 
    plotWidth = 960, plotHeight = 540, progressBar = TRUE) {
  dir.create(path = directory, recursive = TRUE, showWarnings = FALSE)
  timeSeriesCount <- length(timeSeriesList)
  computingCluster <- makeCluster(detectCores())
  registerDoParallel(computingCluster)
  foreach(i = 1:timeSeriesCount, .packages = c("data.table", "ggplot2", "tcltk"),
          .export = "geom_dispatch") %dopar% {
    if (progressBar) {
      if (!exists("progBar")) {
        progBar <- tkProgressBar(max = timeSeriesCount, title = "Plot progress")
      }
      setTkProgressBar(progBar, value = i, label = paste0("Creating plot ", i, "/", timeSeriesCount))
    }
    if (is.null(names(timeSeriesList[[i]]))) { # if multivariate: row names (1st dim)
      if (is.matrix(timeSeriesList[[i]])) {
        timestamps <- 1:nrow(timeSeriesList[[i]])
      } else {
        timestamps <- 1:length(timeSeriesList[[i]])
      }
    } else {
      timestamps <- as.POSIXct(names(timeSeriesList[[i]]), tz = "GMT")
    }
    png(filename = paste0(directory, names(timeSeriesList)[i], ".png"),
        width = plotWidth, height = plotHeight)
    if (is.matrix(timeSeriesList[[i]])) {
      dataTable <- data.table(Time = timestamps, timeSeriesList[[i]])
      dataTable <- melt(dataTable, id.vars = "Time", variable.name = "Attribute",
                        value.name = "Value")
      print(ggplot() +
              geom_dispatch(data = dataTable,
                            mapping = aes(x = Time, y = Value, color = Attribute)) +
              ggtitle(names(timeSeriesList)[i]))
    } else {
      print(ggplot() +
              geom_dispatch(data = data.table(Time = timestamps, Value = timeSeriesList[[i]]),
                            mapping = aes(x = Time, y = Value)) +
              ggtitle(names(timeSeriesList)[i]))
    }
    dev.off()
    return(TRUE)
  }
  stopCluster(computingCluster)
  return(TRUE)
}
