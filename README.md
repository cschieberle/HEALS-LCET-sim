<p align="right">
<img src="http://www.heals-eu.eu/wp-content/uploads/2013/10/logo_heals_spacing.png">
</p>

# Lifecourse Exposure Trajectory Simulation
Lifecourse exposure trajectory simulation for HEALS

The simulation of life-long exposure trajectories is based on the following R library: 
```{r}
library(lifeCourseExposureTrajectories)
```

Please refer to the following source when using the library:

> **Schieberle, C.; Li, N.; Friedrich, R.: *Life course exposure trajectories: A probabilistic model for prospective and retrospective exposure characterization.* 19th International Symposium on Environmental Pollution and its Impact on Life in the Mediterranean, Rome, Italy, 4-6 October 2017.**

## Configuration

To run simulations a configuration object needs to be created:
```{r}
config = lifeCourseExposureTrajectories::defaultConfig()
```

Check the help page on the above command for a list of parameters.
You can always override the (default) configuration by setting configuration parameters individually as follows.

```{r}
# Set the number of master cores which defaults to 2 and should only be raised when you have
# 8 or more cores on your system.
# Make sure that your system has at least 2 additional cores for running exposure sampling.
config[["NUM_MASTER_CORES"]] <- 2

# It is always a good idea to leave at least one core left unused for OS and other issues.
config[["NUM_CORES"]] <- parallel::detectCores() - config[["NUM_MASTER_CORES"]] - 1

# This is used by parallel::makeCluster to set-up a cluster. Default is "PSOCK".
# See help page for details.
config[["CLUSTER_TYPE"]] <- "PSOCK"

# The run-time and memory footprint is largely determined by the number of trajectories to
# be calculated (num.sim) and the number of samples taken per trajectory (sample.size).
config[["NUM_SIM"]] <- num.sim          # e.g. 100
config[["SAMPLE_SIZE"]] <- sample.size  # e.g. 1000

# Shall output files be generated. If set to TRUE, make sure that the output folder is on
# a local drive as opposed to a network share for performance reasons.
config[["WRITE_OUTPUT"]] <- TRUE

# The following path variables should not be set individually but should be defined using the
# parameters of lifeCourseExposureTrajectories::defaultConfig(). See help page for details.
config[["PATH_OUTPUT"]] <- paste0(path, '/', subfolder.output)
config[["PATH_EXPOSURE"]] <- paste0(path, '/', subfolder.exposure)
config[["EMPLOYMENT_MAPPING"]] <- paste0(config[["PATH_EXPOSURE"]], '/', employment.mapping)
```

## General data struture

The expected folder structure is as follows:

* data
  + PM25
  + NO2
  + Cadmium
  + ...
* output-2018-06-05

Within each stressor folder (e.g. PM25) a number of CSV files is expected. Each of them having the following name format:

> PM25_exposure_sample_f4331967-15d1-4094-be76-a34fa77fbf4c.csv

where PM25 is a stressor and f4331967-15d1-4094-be76-a34fa77fbf4c is an identifier of an individual.
