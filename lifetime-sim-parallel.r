options(stringsAsFactors = TRUE)
library(bit64)
library(tidyr)

#detach("package:lifeCourseExposureTrajectories", unload=TRUE)
library(devtools)
devtools::document("..//lifeCourseExposureTrajectories")
devtools::install("..//lifeCourseExposureTrajectories")
library(lifeCourseExposureTrajectories)


FOLDER_CROSS_SECTIONAL <- "K:/EU-SILC cross-sectional/"
FOLDER_LONGITUDINAL <- "K:/EU-SILC longitudinal/"


########################################################################################################
# Step 1: Read cross-sectional and longitudinal EU-SILC data and create the sequence regression tree.
########################################################################################################

# Read cross-sectional data of individuals from H-files.
# First, read household data:
H.files.cs <- paste0(
  FOLDER_CROSS_SECTIONAL,  
  list.files(FOLDER_CROSS_SECTIONAL, pattern="UDB_c(.+)H_ver(.+).csv", recursive = T)
)

cs.hh.data.all <- lifeCourseExposureTrajectories::readHouseholdData(H.files.cs)

# Second, read individual data (P- and R-files) and merge it with household data:
P.files.cs <- paste0(
  FOLDER_CROSS_SECTIONAL, 
  list.files(FOLDER_CROSS_SECTIONAL, pattern="UDB_c(09|10|11|12|13)P_ver(.+).csv", recursive = T)
)
R.files.cs <- paste0(
  FOLDER_CROSS_SECTIONAL, 
  list.files(FOLDER_CROSS_SECTIONAL, pattern="UDB_c(09|10|11|12|13)R_ver(.+).csv", recursive = T)
)

cs.data.all <- lifeCourseExposureTrajectories::readIndividualData(P.files.cs, R.files.cs, cs.hh.data.all)

# Free up household data to reduce memory footprint.
rm(cs.hh.data.all)

# Show countries covered by cross-sectional data:
unique(cs.data.all$COUNTRY)


# Now, read longitudinal data derived from periodic interviews over 4 years.
# Again, read household data first (H-files):
#
H.files.lon <- paste0(FOLDER_LONGITUDINAL, 
                     list.files(FOLDER_LONGITUDINAL, pattern="UDB_(l|L)(06|07|08|09|10|11|12|13)H_ver(.+).csv", recursive = T))

lon.hh.data.all <- lifeCourseExposureTrajectories::readHouseholdData(H.files.lon)

# If only perople from a specific subset of countries are considered, think about
# filtering by country code as this will reduce overall memory usage!
#
# IMPORTANT:  In the following, we will drop all countries but Austria.
#
lon.hh.data.all <- subset(lon.hh.data.all, COUNTRY %in% c("AT"))

# Secondly, all longitudinal data on individuals is collected (P-files):
P.files.lon <- paste0(FOLDER_LONGITUDINAL, 
                      list.files(FOLDER_LONGITUDINAL, pattern="UDB_(l|L)(06|07|08|09|10|11|12|13)P_ver(.+).csv", recursive = T))

data.rshp <- lifeCourseExposureTrajectories::readIndividualSeq(P.files.lon, lon.hh.data.all)
unique(data.rshp$CNTRY)

# Free up household data to reduce memory footprint.
rm(lon.hh.data.all)

#data.seq <- seqdef(data.rshp, 8:11, alphabet = data.alphabet[1:11], states = data.scodes[1:11], labels = data.labels[1:11], xtstep = 1)

# markov.chain <- lifeExposureTrajectories::asMarkovChain(data.seq, data.scodes)
# lifeExposureTrajectories::plotMarkovChain(markov.chain)

data.rshp.subset <- lifeCourseExposureTrajectories::stratifiedSeqSample(data.rshp, 0.25)
data.seq <- seqdef(data.rshp.subset, 8:11, alphabet = data.alphabet[1:11], states = data.scodes[1:11], labels = data.labels[1:11], xtstep = 1)
rm(data.rshp)
gc()

st <- lifeCourseExposureTrajectories::determineSeqTree(data.rshp.subset, dist.method = "OM", dist.indel = 1)

# Install GraphViz if you'd like to see a visualisation of the tree: http://www.graphviz.org
seqtreedisplay(
  st, 
  type = "d", 
  border = NA, 
  only.leaf = F, 
  with.legend = F, 
  gvpath = 'C:\\Program Files (x86)\\GraphViz'
)

edu.data <- lifeCourseExposureTrajectories::getCrossSectEdu(cs.data.all, country.filter = unique(data.rshp.subset$CNTRY))
rm(cs.data.all)
gc()


########################################################################################################
# Step 2: Run the lifecourse trajectory model on a number of individuals in parallel.
########################################################################################################
library(readxl)
library(parallel)
library(tictoc)

# Create a configuration object of the life-course trajectory model and
# set the path to the input files.
#
config <- lifeCourseExposureTrajectories::defaultConfig(
  path = "Y:/Users/xnl/KUNO_kids/woman",
  subfolder.output = paste0("output-", format(Sys.Date(), format="%Y-%m-%d")),
  write.output = TRUE,
  subfolder.exposure = "data",
  employment.mapping = "employment_mapping.xlsx",
  sample.size = 100,
  num.sim = 100
)
config[["CSV_SEP"]] <- ';'
config[["CSV_DEC"]] <- ','
config[["CLUSTER_OUTFILE"]] <- paste0(config[["PATH_OUTPUT"]], "/clusterout.txt")

# The list of stressors for which the model will run.
#
#config[["stressors"]] <- c("NO2", "UV", "EMF")
config[["stressors"]] <- c("PM25")

# Read information of the individuals for which the lifecourse exposure
# should be modelled.
#
individuals_csv <- read.csv(
  paste0("N:/tfu/552_HEALS/Projektarbeit/Application/KUNO_Kids/KUNO_trial/KUNO_women_SES_variables.csv"),
  sep = ";"
)
individuals <- NULL
for (i in 1:nrow(individuals_csv)) {
  if (individuals_csv[i,]$edcat_kuno > -1) {
    individuals <- c(
      individuals,
      new("Individual",
          id = as.character(individuals_csv[i,]$identifier),
          age = individuals_csv[i,]$age,
          sex = "F",
          edulevel = (individuals_csv[i,]$edcat_kuno - 1)
      )
    )
  }
}

# Create a cluster object. Here a local cluster with 2 master cores will be initialized.
# The remaining cores will be used to run the exposure estimation and sampling.
#
cl.master <- makeCluster(config[["NUM_MASTER_CORES"]], outfile = config[["CLUSTER_OUTFILE"]])

# Make sure all nodes in the cluster hold the necessary R libraries.
# Furthermore, share some variables across the cluster.
#
clusterEvalQ(cl.master, library(lifeCourseExposureTrajectories, quietly = T))
clusterEvalQ(cl.master, library(readxl, quietly = T))
clusterEvalQ(cl.master, library(parallel, quietly = T))
clusterExport(cl.master, varlist = c("config", "st", "data.seq", "edu.data", "individuals"))

# Increase memory limit ...
#
memory.limit(24000)

# Start a timer ...
#
tic()

# Run the estimation and parallelize across the individuals.
#
sim.results <- parLapply(
  cl.master,
  1:10, 
  par_lifeCourseExposure
)

# Combine the statistics of all individuals.
#
sample.exp.stats <- do.call(rbind.data.frame, sim.results)

# Stop the timer and shut down the cluster.
#
exectime <- toc()
exectime <- exectime$toc - exectime$tic
exectime

stopCluster(cl.master)  

length(unique(sample.exp.stats$INDIV_SUBJID))
