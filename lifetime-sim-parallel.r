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

################################
#
# PARALLEL PART
#
################################
library(readxl)
library(parallel)

config <- lifeCourseExposureTrajectories::defaultConfig(
  path = "N:\\tfu\\552_HEALS\\Projektarbeit\\WP11",
  sample.size = 100
)

config[["NUM_SIM"]] <- 100
config[["stressors"]] <- c("NO2", "UV", "EMF")

individuals <- read.csv(paste0(config[["PATH"]], "\\Stream 5 data\\stream5SampleData3.csv"))
stopifnot( length(individuals$id) == length(unique(individuals$id)) )


cl.master <- makeCluster(config[["NUM_MASTER_CORES"]])

clusterEvalQ(cl.master, library(lifeCourseExposureTrajectories))
clusterEvalQ(cl.master, library(readxl))
clusterEvalQ(cl.master, library(parallel))

clusterExport(cl.master, "config")
clusterExport(cl.master, "st")
clusterExport(cl.master, "data.seq")
clusterExport(cl.master, "edu.data")
clusterExport(cl.master, "individuals")

require(tictoc)
tic()

sim.results <- parLapply(
  cl.master,
  402:403, 
  par_lifeCourseExposure
)
sample.exp.stats <- do.call(rbind.data.frame, sim.results)

exectime <- toc()
exectime <- exectime$toc - exectime$tic
exectime

stopCluster(cl.master)  

unique(sample.exp.stats$INDIV_SUBJID)



##################
##################
# 
# names(d)
# 
# sort(unique(sample.exp.all$age))
# 
# # name critical life stages ...
# sample.exp.all$CRITICAL_LIFE_STAGE <- ""
# sample.exp.all[ sample.exp.all$age <= 3, ]$CRITICAL_LIFE_STAGE <- "Infancy"
# sample.exp.all[ sample.exp.all$age > 3 & sample.exp.all$age <= 11, ]$CRITICAL_LIFE_STAGE <- "Childhood"
# sample.exp.all[ sample.exp.all$age > 11 & sample.exp.all$age <= 17, ]$CRITICAL_LIFE_STAGE <- "Adolescence"
# sample.exp.all[ sample.exp.all$age > 17 & sample.exp.all$age <= 39, ]$CRITICAL_LIFE_STAGE <- "Adulthood.below.40"
# sample.exp.all[ sample.exp.all$age > 39 & sample.exp.all$age <= 64, ]$CRITICAL_LIFE_STAGE <- "Adulthood.below.65"
# sample.exp.all[ sample.exp.all$age > 64, ]$CRITICAL_LIFE_STAGE <- "Adulthood.65.and.older"
# 
# # ... and generate aggregate stats
# sample.exp.stats <- aggregate(value ~ INDIV_SUBJID + STRESSOR + CRITICAL_LIFE_STAGE, data=sample.exp.all, FUN=quantile, probs=c(0.025, 0.5, 0.975))
# sample.exp.stats <- cbind(sample.exp.stats[, 1:3 ], data.frame(sample.exp.stats$value))
# names(sample.exp.stats) <- c("INDIV_SUBJID", "STRESSOR", "CRITICAL_LIFE_STAGE", "value.2.5PCT", "value.50PCT", "value.97.5PCT")
# 
# sample.exp.stats <-
# merge(
#   merge(
#     merge(
#       sample.exp.stats[ sample.exp.stats$CRITICAL_LIFE_STAGE == "Infancy", ],
#       sample.exp.stats[ sample.exp.stats$CRITICAL_LIFE_STAGE == "Childhood", ],
#       by = c("INDIV_SUBJID", "STRESSOR"),
#       suffixes = c(".Infancy", ".Childhood")
#     ),
#     merge(
#       sample.exp.stats[ sample.exp.stats$CRITICAL_LIFE_STAGE == "Adolescence", ],
#       sample.exp.stats[ sample.exp.stats$CRITICAL_LIFE_STAGE == "Adulthood.below.40", ],
#       by = c("INDIV_SUBJID", "STRESSOR"),
#       suffixes = c(".Adolescence", ".Adulthood.below.40")
#     ),
#     by = c("INDIV_SUBJID", "STRESSOR"),
#     suffixes = c("", "")
#   ),
# 
#   merge(
#     sample.exp.stats[ sample.exp.stats$CRITICAL_LIFE_STAGE == "Adulthood.below.65", ],
#     sample.exp.stats[ sample.exp.stats$CRITICAL_LIFE_STAGE == "Adulthood.65.and.older", ],
#     by = c("INDIV_SUBJID", "STRESSOR"),
#     suffixes = c(".Adulthood.below.65", ".Adulthood.65.and.older")
#   ),
#   by = c("INDIV_SUBJID", "STRESSOR"),
#   suffixes = c("", "")
# )
# 
# 
# sample.exp.stats.plot.ypos <- aggregate(EXP_97.5PCT ~ INDIV_SUBJID + STRESSOR, data=d, FUN=max)
# names(sample.exp.stats.plot.ypos) <- c("INDIV_SUBJID", "STRESSOR", "PLOT.YPOS")
# sample.exp.stats.plot.ypos$PLOT.YPOS <- 0.9 * sample.exp.stats.plot.ypos$PLOT.YPOS
# 
# sample.exp.stats <- merge(sample.exp.stats, sample.exp.stats.plot.ypos, by=c("INDIV_SUBJID", "STRESSOR"))
# 
# sort(unique(d$STRESSOR))
# 
# 
# library(ggplot2)
# 
# stressor_unit_labeller <- c(
#   PM25 = "PM2.5\n[ug.d-1]",
#   NO2 = "NO2\n[ug.d-1]",
#   SOMO35 = "SOMO35\n[ug.y-1]",
#   UV = "UV\n[J.cm-2.d-1]",
#   EMF = "EMF\n[mG.d-1]",
#   Chromium = "Chromium\n[ug.kgbw-1.d-1]",
#   Lead = "Lead\n[ug.kgbw-1.d-1]",
#   Mercury = "Mercury\n[ug.kgbw-1.d-1]",
#   Mould = "Mould\n[unitless]",
#   DEHP = "DEHP\n[ng.d-1]",
#   DINP = "DINP\n[ng.d-1]",
#   DIDP = "DIDP\n[ng.d-1]"
# )
# 
# #
# 
# ggplot(
#   #data = merge(d, d.plot.theme, by="STRESSOR"),
#   #data = merge(merge(d, d.plot.theme, by="STRESSOR"), sample.exp.stats, by=c("INDIV_SUBJID", "STRESSOR")),
#   data = merge(d, sample.exp.stats, by=c("INDIV_SUBJID", "STRESSOR")),
#   aes(x = AGE, y = EXP_MEDIAN, ymin = EXP_25PCT, ymax = EXP_75PCT)) + 
#   geom_ribbon(alpha = 0.3, aes(fill = STRESSOR)) + 
#   geom_line(size = rel(1.2), aes(colour = STRESSOR)) + 
#   geom_line(aes(x = AGE, y = EXP_97.5PCT, colour = STRESSOR), linetype = "longdash") +
#   geom_line(aes(x = AGE, y = EXP_2.5PCT, colour = STRESSOR), linetype = "longdash") +
#   geom_vline(xintercept = 3, colour="black", linetype = "solid") +
#   geom_vline(xintercept = 12, colour="black", linetype = "solid") +
#   geom_vline(xintercept = 18, colour="black", linetype = "solid") +
#   geom_vline(xintercept = 40, colour="black", linetype = "solid") +
#   geom_vline(xintercept = 65, colour="black", linetype = "solid") +
#   theme_bw() + 
#   # geom_label(aes(x=1.5,  y=PLOT.YPOS, label=paste0(
#   #   signif(value.50PCT.Infancy, 4), " [",
#   #   signif(value.2.5PCT.Infancy, 4), "; ",
#   #   signif(value.97.5PCT.Infancy, 4), "]"
#   # ))) +
#   # geom_label(aes(x=7.5,  y=PLOT.YPOS, label=paste0(
#   #   signif(value.50PCT.Childhood, 4), " [",
#   #   signif(value.2.5PCT.Childhood, 4), "; ",
#   #   signif(value.97.5PCT.Childhood, 4), "]"
#   # ))) +
#   # geom_label(aes(x=15,  y=PLOT.YPOS, label=paste0(
#   #   signif(value.50PCT.Adolescence, 4), " [",
#   #   signif(value.2.5PCT.Adolescence, 4), "; ",
#   #   signif(value.97.5PCT.Adolescence, 4), "]"
#   # ))) +
#   # geom_label(aes(x=24.5,  y=PLOT.YPOS, label=paste0(
#   #   signif(value.50PCT.Adulthood, 4), " [",
#   #   signif(value.2.5PCT.Adulthood, 4), "; ",
#   #   signif(value.97.5PCT.Adulthood, 4), "]"
#   # ))) +
#   theme(legend.position="none") +
#   labs(title = paste0("exposure of ", INDIV_SUBJID), x = "age in years", y = "") +
#   #scale_x_continuous(breaks = sort(unique(d$AGE))) +
#   scale_x_continuous(breaks = c(3,12,18,40,65)) +
#   facet_grid(STRESSOR ~ ., scales="free", labeller = labeller(STRESSOR = stressor_unit_labeller))
# 
# # ####
# # 
# # sim.result.seq <- seqdef(
# #   sim.results.rshp, 
# #   2:82, 
# #   alphabet = data.alphabet, 
# #   states = data.scodes, 
# #   labels = data.labels, 
# #   xtstep = 1,
# #   cpal = c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(length(data.scodes) - 12, "Set1"))
# # )
# # 
# # dev.off()
# # seqdplot(sim.result.seq, withlegend = FALSE, border = NA)
# # #seqIplot(sim.result.seq, withlegend = FALSE, border = NA)

