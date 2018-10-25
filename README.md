# ImpactAggregationClustering

[![R 3.4.3](https://img.shields.io/badge/R-3.4.3-blue.svg)](https://www.r-project.org/)
[![packrat 0.4.8.1](https://img.shields.io/badge/packrat-0.4.8.1-blue.svg)](https://rstudio.github.io/packrat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the code to reproduce the experiments and analyses of the paper

> Holger Trittenbach, Jakob Bach, Klemens Böhm, "Understanding the Effects of Energy-Data Aggregation on Clustering Quality"

For more information on our research project, see the companion website: https://www.ipd.kit.edu/clustagg/.

## Setup

The implementation requires [R >= 3.4.3](https://cran.r-project.org/).
Dependencies are managed with [packrat](https://rstudio.github.io/packrat/).
To download the required packages, run the following command inside the project directory:

```bash
$ Rscript -e 'packrat::restore()'
```

The package `FastTSDistances` must be installed manually from a separate [repository](https://github.com/Jakob-Bach/FastTSDistances).

## Overview

The following table provides a short overview of the R files and their purpose.
The files in sub-directories are explained in separate sections.

The main experiment script is located at `batchscripts/SparkAggregationBasedScript.R`.
Files ending with `Utility` only contain functions or global variables; they can be sourced without triggering heavy computations.

|File|Content|
|---|---|
|.Rprofile|Calls a *packrat* initialization routine when you open our project.|
|AggregationUtility.R|Contains methods for uni- and multi-variate (e.g., min and max for each time interval) piecewise aggregation of uni- and multi-variate sequences. There is a default list of aggregations, which you could extend for your experiments. There also is a [batchtools](https://mllg.github.io/batchtools/articles/batchtools.html) routine to aggregate one or multiple sequences.|
|BatchExperimentEvaluationUtility.R|Contains evaluation functions for the results of a batch experiment. This includes generating plots for single clusterings and forecasts as well as creating summary tables.|
|`batchscripts/`|Contains scripts which run without manual intervention, e.g., the main experiment, plotting or reporting. Necessary configuration options can be passed as parameters.|
|BatchUtility.R|Contains functions for the different steps in a batch experiment. This includes parameter checking, pre-processing, aggregation, clustering and evaluation (including forecasting). This file mostly contains high-level functionality and calls methods from utility files.|
|ClusteringEvaluationUtility.R|Contains methods to evaluate a clustering result (cluster assignments, dissimilarity matrices). Possible ways to evaluate are 1) external cluster validity indices comparing results on different aggregation levels, 2) external validity indices comparing against a pre-defined ground truth vector (like machines, days of week etc.), 3) internal validity indices, 4) cluster size entropy, 5) cluster counts, 6) forecasting on the raw datasets using clustering results and 7) forecasting on aggregated datasets using clustering results.|
|ClusteringUtility.R|Contains a [batchtools](https://mllg.github.io/batchtools/articles/batchtools.html) routine for clustering (considering multiple datasets, aggregation functions, aggregation levels, dissimilarities, clustering algorithms) and two methods for cluster assignment. Clustering algorithms are wrapped into methods in a standardized way. Thus, you could also extend the default list for your own experiments. There also are functions to create cluster plots and to copy plots of sequences according to their cluster membership.|
|DissimilarityUtility.R|Contains lists of uni- and multivariate dissimilarities which are used by default in the experiments. Similar to the lists of aggregations and clustering algorithms, you could also add your own dissimilarity functions.|
|GeneralPlotRoutinesUtility.R| Contains methods to plot uni- as well as multi-variate sequences.|
|PreprocessingUtility.R|Contains functions to pre-process datasets. This includes processing of CSV files, interpolating, conversion of tables to sequence lists and summarizing multiple data files from the same directory (missing values, constant series etc.).|
|`packrat/`|Contains *packrat* configuration files and also stores the downloaded dependencies after initializing *packrat*.|
|`preprocscripts/`|Contains scripts for additional pre-processing, e.g., outlier removal.|
|`presentationscripts/`|Contains scripts to analyze data and reproduce plots from the paper.|

## Batch Scripts (directory : `batchscripts`)

All files in the `batchscripts` directory are non-interactive, and log progress information to the console.
Depending on the script, named or unnamed parameters must be passed upon execution.
To avoid path conflicts, execute these scripts from the main project directory, not the `batchscripts` sub-directory.
Example:

```
path/to/this/project>Rscript batchscripts/<<Script name>>.R <<params>>
```

### `SparkAggregationBasedScript.R`

This script is the main analysis pipeline from pre-processing to clustering to evaluation.
Although we used pre-aggregated data (created by *Spark*) in our experiments, the script can also execute the aggregation locally.
The input directory for an experiment is passed as the `dataPath` parameter.
This directory must contain textfiles in CSV format with the columns *machine_name*, *sensor_name*, *date_start*, *date_end*, and further columns explained below.
_date_ columns expect timestamps in the format `2017-04-01T00:00:00.000+02:00`.

#### Pre-Processing

If the `dataPath` contains only one CSV file, aggregations are created locally at predefined levels (currently hard-coded: 30 sec sampling interval, 1 min, 5 min, 10 min, 15 min, 30 min, 1 h, 2 h, 6 h).
**Additional Columns:** The file must contain an additional column *mean* which contains values to be aggregated.

One can provide multiple CSV files where each contains data for a specific aggregation level.
File names must contain exactly one substring in the format `[0-9]+_[a-z]+`, e.g., `15_min`.
This is only to distinguish the files; the actual sampling rate is inferred automatically from the data.
**Additional Columns:** Besides the standard identifier columns, the input file must have a column *mean* and several additional (arbitrarily named) columns with aggregated values, e.g., *median*.
The *mean* column is used for summary operations before creating sequences, e.g., detecting constant sequences, and is considered the base-level data.

The sequences for each aggregation function and level are saved in `aggData/`.
As a special case, the *mean* is only saved for the largest dataset file (which is assumed to have the highest sampling rate) and considered the base-level data (see the paper for details on the base-level).
There also are overview tables that describe the sequences: `machineSensorDayCombinations.rds` (fixed-length sequences) and `timeSeriesOverviewTable.rds` (variable-length sequences).
See the parameter `intervalLength` for details on how to use these two options.
The dissimilarities and clustering algorithms in the analysis cannot be configured via console parameters.
To add dissimilarities or clustering algorithms, update the default algorithm lists in `ClusteringUtility.R` and `DissimilarityUtility.R`.

#### How to extended the pre-processing

Custom pre-processing steps, e.g., outlier removal, require an additional script which
1. execute the first part of the existing pre-processing pipeline (read CSVs, serialize as `.rds` files).
2. executes the custom pre-processing.
3. saves the resulting *data.table* in an RDS file.

The pipeline of `SparkAggregationBasedScript.R` can pick up intermediate results and directly start with the second part of pre-processing (conversion to sequences).
See `preprocscripts/` for examples.

#### Details on the analysis pipeline

The script `SparkAggregationBasedScript.R` executes the following steps:
* Save experiment parameters to `params.rds` and `params.txt` under `results/`.
* Read input CSVs, save them as RDS files with some additional columns and simple data type conversions. This is skipped if already done; for our experiments we have provided the RDS directly after additional pre-processing, see `preprocscripts/`.
* Convert input *data.table*s to lists of sequences and save them to `aggData/` (or four different sub-directories if forecasting enabled). This is skipped if `aggData` already contains the aggregated sequences.
* Perform local aggregation if only one CSV is provided.
* Run the clustering: Save dissimilarity matrices to `dissMatrices/` and the cluster assignments to `results/`. This step requires the most execution time.
* Run the standard evaluation: Results are saved as *data.table*s to `results/`.
  * cluster count (not described in the paper, but part of our published experiments)
  * cluster size entropy (not described in the paper, but part of our published experiments)
  * internal validity
  * external validity
  * ground truth validity: currently supported ground truth assignments are machines, sensors, days, time of days, day of week, weekdays vs weekend (not described in our paper, but part of our published experiments)
* Run the forecasting evaluation (optional): Assign test data to training data clusters (time consuming, due to dissimilarity computations) and make instance-based forecasts. Save outputs to `testDissMatrices/` and `testDataAssignmentResult.rds` in `results/`. Compute base and relative forecast error and save them as *data.table*s. (not described in our paper, not part of our published experiments)
* Execute `BatchExperimentReportScript.R` and `BatchExperimentPlotScript.R` to create summary tables and plots.

#### Script parameters

The parameters of `SparkAggregationBasedScript.R` are:

- `dataPath`: Absolute or relative path to the directory containing the input data files (CSV or already RDS). The results of the analysis will be saved in sub-directories. Default:  `./`. If your path contains `data`, as in `.../data/MyDataset/`, standard evaluation plots created by internally calling `BatchExperimentPlotScript.R` will be saved to `.../plots/MyDataset/`, else they also go into the `dataPath` directory.
- `plotPath`: The path to create cluster plots. Each cluster can be represented by one plot showing all its members and the centroid (if existing). However, considering the cross-product of different experiments parameters (aggregations, aggregation levels, dissimilarities, clustering algorithms), using the plot functionality a) takes a lot of time and b) consumes a lot of space on your hard drive. It is rather recommended to create plots for some (interesting) experiment settings afterwards. You can disable plotting by setting the value to `""`, which is also the default. Note that this does not affect the few evaluation plots created by `BatchExperimentPlotScript.R`.
- `datasetName`: String which will be used to name files storing the lists of sequences, dissimilarity matrix files and the dataset in all result tables. Because of internal naming conventions, this parameter must not contain an underscore. Default: `MyDataset`.
- `normalize`: Should all sequences be normalized? Supported are "zscore", "minmax" (both applied to each sequence individually) or any other value, which results in no normalization (default).
- `testRatio`: Allows enabling and disabling forecast evaluation. If `0` (default), forecasting evaluation is switched off, else the dataset it split (chronologically) with a fraction of `testRatio` belonging to the test set. Instead of one directory `aggData`, four folders are used to save the lists of sequences: On the one hand, there are training and test set, one the other hand, you have to distinguish between input data set (relevant for clustering, assignment) and full dataset (relevant for forecasting, contains the input sequences on which the forecasts are made as well as the target sequences which should be forecast). Forecasting is not supported for sequences of variable length.
- `aggregations`: Aggregation functions to be used in the analysis. Have to correspond to columns in the provided input dataset files. The aggregation string consists of different aggregations separated by the pipe character. You can also use multi-variate aggregations denoted by an underscore like `min_max`. Default: `mean|median|min|max|stddev|skewness|kurtosis`. If you want to perform local aggregation, you can use these seven names; multi-variate aggregations as combinations of them are possible as well.
- `machines`: Use a subset of the machines contained in the dataset. Different machines are separated by the pipe character, e.g., `AVT_01|AVT_04|AVT_06`. Default: `""` (no sub-setting).
- `intervalLength`: Our analysis supports sequences of fixed as well as variable length.
  - Positive values of this parameter indicate the length of sequences in seconds (not time steps/elements!). The default is `86400`, i.e., one day, but you can also use fractions of this (like 6 hours), as long as they divide one day without remainder. Sequences are created for each machine and sensor, automatically removing sequences where the sensor only has one value (constant sequences) or where 1% or more of the sequence are missing.
  - The second alternative is creating sequences of variable length. Variable length sequences start if the value of a sensor switches from a base value (like 0) to another value and they end if the sensor comes back to the base value. However, a sensor might also switch to the base value only for a very short time during a period of activity. If you pass a negative value for `intervalLength`, this will be interpreted as the minimum time a sensor needs to remain on the base value (or not send any data at all) to separate two sequences. A minimum sequence length of 30 min is enforced in the code (hard-coded). All values which make up 5% or more of a sensor's measurements qualify as base value. If multiple values fulfill this criterion, a warning is printed and the whole range from the lowest to the highest base value is considered as base (reasonable for some sensors). If this does not yield sensible results, consider the third possibility.
  - As a third possibility, you can pass an `intervalLength` of `0`. This denotes the case that you want variable length sequences, but use the sequence starts and ends from another dataset (or defined manually). In this case, you need to place a file called `timeSeriesOverviewTableInput.rds` in the directory `results/` under the `dataPath`. Such a table could for example be obtained from a different experiment with variable length sequences (sequence definitions are saved in a file called `timeSeriesOverviewTable.rds`). Such a table needs to contain the columns
    - `machine_name`, `sensor_name` (character/factor)
    - `start_time_numeric`, `end_time_numeric` (numeric)
    - `startDay` (character like `2017-12-31`)
    - `startTimeOfDay` (integer, denoting the number of the sequence on this particular machine, sensor, day)

All parameters are passed with their name, the order does not matter.
All parameters are optional (as they all have defaults).
The script checks if the parameters' values are reasonable and might change them to defaults or even stop if they are not.

Example calls:

```bash
>Rscript batchscripts/SparkAggregationBasedScript.R --dataPath="../data/MyExperiment/" --datasetName="ExpData" --intervalLength=-1800 --machines="AVT_01|AVT_04|AVT_06"

>Rscript batchscripts/SparkAggregationBasedScript.R --dataPath="../data/MyExperiment/" --datasetName="ExpData" --intervalLength=43200 --testRatio=0.3 --normalize --aggregations="mean|mean_min_max|median_skewness"
```

### `BatchExperimentPlotScript.R`

This script creates plots based on the results of a batch experiment.
For several evaluation categories, it plots all values of a selected evaluation index (like silhouette as an internal validity index) against aggregation levels, together with box plots.
Furthermore, it plots an aggregate (currently implemented: the mean) of these evaluation measures grouped by a category, e.g., to compare the mean base level internal validity by dissimilarity and aggregation level.
Evaluation categories which are currently considered are:

- cluster count
- cluster size entropy
- internal validity indices (base, current), currently using silhouette as index
- external validity indices comparing aggregation levels (base, previous), currently using inverted van Dongen as index
- all kinds of ground truth evaluation (from theoretical point of view also external validity), currently using inverted van Dongen as index
- (if forecasting evaluation used) relative and base forecast errors for the weighted clustering forecast, expressed as ratio against the naive forecast

The script has only one parameter which is passed unnamed and corresponds to the `dataPath` in the main experiment script.
It creates a `plotPath` from that by replacing `data` in the path with `plots`.

Example call:

```bash
>Rscript batchscripts/BatchExperimentPlotScript.R ../data/MyExperiment/
```

### `BatchExperimentReportScript.R`

This script creates a report based on the results of a batch experiment.
It prints the standard R *summary()*, checks for NAs and creates tables comparing the mean of a certain evaluation index for different experiment categories (like clustering, dissimilarity etc.).
Currently considered in the report are:

- cluster count
- cluster size entropy
- internal validity indices (base, current), currently using silhouette as index
- external validity indices comparing aggregation levels (base, previous), currently using inverted van Dongen as index
- all kinds of ground truth evaluation (from theoretical point of view also external validity), currently using inverted van Dongen as index
- relative and base forecast errors expressed as ratio against the naive forecast

It writes an `evaluation.txt` file to the `results/` directory of your experiment, using *Markdown*-style headings (but tables are included as produced by R).
It only has one parameter which is passed unnamed and corresponds to the `dataPath` in the main experiment script.
Although the file works as batch script, you can also use it interactively and run only selected parts, e.g., if you want to drill down into the results or create cluster/forecast plots (calls to the corresponding routines are included in the script, but commented).

Example call:

```bash
>Rscript batchscripts/BatchExperimentReportScript.R ../data/MyExperiment/
```

### `ForecastingScript.R`

This script only executes the steps belonging to forecasting evaluation, i.e., cluster assignments, base forecasting and relative forecasting.
It might be called if an error occurs during assignment or the assignment strategy has changed.
The script assumes the same pre-processing and clustering as performed in `SparkAggregationBasedScript.R`.
In fact, it mainly copies the second part of that script (leaving out unnecessary parts, reading existing results from previous steps as RDS files, even the input parameters).
It only has one parameter which is passed unnamed and corresponds to the `dataPath` in the main experiment script.
When running this script, a hint is added to the `params.txt` file of the original experiment, stating that the forecasting part has been re-executed.

Example call:

```bash
>Rscript batchscripts/ForecastingScript.R ../data/MyExperiment/
```

### `ForecastingTablesScript.R`

This script re-generates the base forecasting and relative forecasting error tables.
It might be called if only the forecasting routine has changed (e.g., a new forecasting technique added), but the rest of the experiment can remain as-is.
The script assumes the same pre-processing, clustering and cluster assignment as performed in `SparkAggregationBasedScript.R`.
In fact, it mainly copies a part of the forecasting evaluation routine (leaving out the assignment phase, reading existing results from previous steps as RDS files).
It only has one parameter which is passed unnamed and corresponds to the `dataPath` in the main experiment script.
When running this script, a hint is added to the `params.txt` file of the original experiment, stating that the forecasting table computation has been re-executed.

Example call:

```bash
>Rscript batchscripts/ForecastingTablesScript.R ../data/MyExperiment/
```

### `MergeSummerWinterCSVScript.R`

The *Spark* aggregation we used for our datasets creates intervals of fixed length regarding a fixed time zone.
When the clock changes to daylight-saving time in March, the time zone of the timestamps in the exported CSV is adapted, e.g., it changes from +01 to +02.
However, the aggregation intervals are still aligned to the original time zone, e.g., if they start at 00:00, 02:00, 04:00 in winter, they start at 01:00, 03:00, 05:00 in summer (so even the interval during which the clock is turned forward spans the same time interval of two hours regarding a "neutral" time zone like UTC).
However, we want to have an alignment which always starts at midnight of the local/relevant time zone.
The current solution is to export two files for aggregates with a sampling interval higher than one hour, one aligned to midnight/"ordinary" time and one aligned to midnight/daylight-saving time.
`MergeSummerWinterCSVScript.R` merges them to one CSV file, taking the part from each file which is correctly aligned.

Example call:

```bash
>Rscript batchscripts/MergeSummerWinterCSVScript.R ../data/MyExperiment/
```

### `RandomDataGenerationScript.R`

This script generates normally distributed random data and computes aggregates on them.
This can be used to analyze if the effects of aggregation can even be found for random data.
The output is one RDS file containing the original data and several RDS files containing the aggregate tables, so they can be used within the analysis pipeline.
The script has the following parameters:

- `dataPath`: Absolute or relative path to the directory where the aggregate table files should be placed. The original table will be saved in a sub-directory. Default:  `./`.
- `days`: Number of days. Default: 60.
- `machines`: Number of machines. Default: 10.
- `samplingInterval`: Interval between measurements in the original data (in seconds). Default: 3.

Arguments are not checked for sanity.
Besides the seven standard aggregates, a column containing a random sample from the aggregation interval of the original data (`sample`) and a column containing random numbers generated independently for each aggregation level (`random`) are added.
If they should be considered in the aggregation script, make sure to include them in the `aggregations` parameter.

Example call:

```bash
>Rscript batchscripts/RandomDataGenerationScript.R --dataPath="../data/MyExperiment/" --machines=5
```

## Pre-Processing Scripts (directory : `preprocscripts`)

Although we have a fully automated analysis pipelines, some datasets might need manual intervention, in particular in the pre-processing step.
The scripts in this directory provide solutions for specific datasets used in our experiments.
They are not suitable for general-purpose use.
Data paths are hard-coded in these scripts.

### `ActivePowerCSVPreprocScript.R`

In the *Active Power* (Jan-Aug 2017) dataset, some machines have a permanently negative value.
This should not happen if amperage as well as voltage are non-negative and simply was a measurement error.
The resulting sequences look like the ones of other machines, only mirrored at the value 0.
To enable comparison of all machines, we add a post-processing step after the usual CSV to RDS conversion:
We invert data for the affected machines (assuming all original measurements have the same sign, we can also determine the aggregates after switching the sign) and save them again as RDS files with the suffix `_valid`.

### `CreateActivePowerJanAugOverviewTable.R`

Script which defines the sequence starts and ends for the *Active Power* (Jan-Aug 2017) dataset with variable length sequences.
We have a simple routine which tries to find the baseline in each dataset, but it did not work for one of our machines, as two baselines were found and one of them was clearly to high.
As a consequence, we define the baseline for this machine manually and merge with the default result obtained for the other machines.
The resulting overview table with the sequence definitions is saved as RDS and can be used with the main experiment script (`--intervalLength=0`).
In fact, we also use it to define variable length sequences for other attributes than power (which do not show a baseline behavior at all, i.e., frequency, positive energy and voltage).

### `CreatePowerFactorJanAugOverviewTable.R`

Script which defines the sequence starts and ends for the *Power Factor* (Jan-Aug 2017) dataset with variable length sequences.
For several machines, two baseline values (0 and 1), which could be used for sequence extraction, are found.
We decide on using 1 as common baseline for all machines.
The resulting overview table with the sequence definitions is saved as RDS and can be used with the main experiment script (`--intervalLength=0`)

### `FrequencyCSVPreprocScript.R`

The *Frequency* (Jan-Aug 2017) dataset generally has values fluctuating around 50 Hz.
However, from time to time, there are sharp drops or the sensor is switched off.
Clustering without additional pre-processing mainly isolates these special sequences.
To obtain further insights, we add a post-processing step after the usual CSV to RDS conversion:
We remove the days which show this unwanted behavior from all aggregation datasets and save the reduced dataset again as RDS files with the suffix `_valid`.

### `MeanAmperageCSVPreprocScript.R`

The *MeanAmperage* (Jan-Aug 2017) dataset generally has a base value behavior, i.e., sensors start and end with a value of zero each day and have positive values in between.
However, sometimes sequences show a different shape, as some machines run over several days.
Apart from that, on some days there are spikes for only very short periods of time, so probably outliers.
To make the dataset more homogeneous regarding base behavior, we add a post-processing step after the usual CSV to RDS conversion:
We remove the days which show this unwanted behavior from all aggregation datasets and save the reduced dataset again as RDS files with the suffix `_valid`.

### `MeanVoltageCSVPreprocScript.R`

The *MeanVoltage* (Jan-Aug 2017) dataset generally has values fluctuating around 400 V.
However, from time to time, there are sharp drops or the sensor is switched off.
Clustering without additional pre-processing mainly isolates these special sequences.
To obtain further insights, we add a post-processing step after the usual CSV to RDS conversion:
We remove the days which show this unwanted behavior from all aggregation datasets and save the reduced dataset again as RDS files with the suffix `_valid`.

### `PositiveEnergyCSVPreprocScript.R`

The *Positive Energy* (Jan-Aug 2017) dataset should only contain increasing values according to a domain expert, as the total energy consumption (in 1 kWh steps) is measured.
However, for some sensors the energy drops to 0 if the sensor is switched off or has drops during periods with otherwise constant behavior.
As a consequence, we add a post-processing step after the usual CSV to RDS conversion:
We remove the days which show decreasing behavior from all aggregation datasets and save the reduced dataset again as RDS files with the suffix `_valid`.

### `PowerFactorCSVPreprocScript.R`

The *PowerFactor* (Jan-Aug 2017) dataset generally has a base value behavior, i.e., sensors start and end with a value of 1 each day and have smaller values (between -1 and 1) during the day.
However, sometimes sequences show a different shape, as some machines run over several days.
Apart from that, on some days there are spikes for only very short periods of time, so probably outliers.
To make the dataset more homogeneous regarding base behavior, we add a post-processing step after the usual CSV to RDS conversion:
We remove the days which show this unwanted behavior from all aggregation datasets and save the reduced dataset again as RDS files with the suffix `_valid`.

## Presentation Scripts (directory : `presentationscripts`)

This directory contains scripts which perform data analysis and create plots for presentation purposes.
The scripts do not run fully automated, but contain hard-coded paths and are partially tied to the datasets from our experiments.

|File|Content|
|---|---|
|EvaluationDemo.R|File which demonstrates evaluation functions on the data from our experiments.|
|PaperPlots.R|File which was used to create plots for our paper.|
|PaperPlotUtility.R|File which provides functions for the `PaperPlots.R.` files.|
