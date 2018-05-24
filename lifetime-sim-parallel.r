options(stringsAsFactors = TRUE)
library(bit64)
library(tidyr)

detach("package:lifeCourseExposureTrajectories", unload=TRUE)
library(devtools)
devtools::document("..//lifeCourseExposureTrajectories")
devtools::install("..//lifeCourseExposureTrajectories")
library(lifeCourseExposureTrajectories)

H.files.cs <- c(
  "K:/EU-SILC cross-sectional/1 C-2004/UDB_c04H_ver 2004-4 from 01-08-09.csv",
  "K:/EU-SILC cross-sectional/2 C-2005/UDB_c05H_ver 2005-5 from 01-08-09.csv",
  "K:/EU-SILC cross-sectional/3 C-2006/UDB_c06H_ver 2006-4 from 01-03-10.csv",
  "K:/EU-SILC cross-sectional/4 C-2007/UDB_c07H_ver 2007-6 from 01-08-11.csv",
  "K:/EU-SILC cross-sectional/5 C-2008/UDB_c08H_ver 2008-7 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/6 C-2009/UDB_c09H_ver 2009-7 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/7 C-2010/UDB_c10H_ver 2010-6 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/8 C-2011/UDB_c11H_ver 2011-5 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/9 C-2012/UDB_c12H_ver 2012-3 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/10 C-2013/UDB_c13H_ver 2013-2 from 01-08-15.csv"
)
cs.hh.data.all <- lifeCourseExposureTrajectories::readHouseholdData(H.files.cs)

P.files.cs <- c(
  "K:/EU-SILC cross-sectional/6 C-2009/UDB_c09P_ver 2009-7 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/7 C-2010/UDB_c10P_ver 2010-6 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/8 C-2011/UDB_c11P_ver 2011-5 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/9 C-2012/UDB_c12P_ver 2012-3 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/10 C-2013/UDB_c13P_ver 2013-2 from 01-08-15.csv"
)
R.files.cs <- c(
  "K:/EU-SILC cross-sectional/6 C-2009/UDB_c09R_ver 2009-7 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/7 C-2010/UDB_c10R_ver 2010-6 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/8 C-2011/UDB_c11R_ver 2011-5 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/9 C-2012/UDB_c12R_ver 2012-3 from 01-03-15.csv",
  "K:/EU-SILC cross-sectional/10 C-2013/UDB_c13R_ver 2013-2 from 01-08-15.csv"
)

cs.data.all <- lifeCourseExposureTrajectories::readIndividualData(P.files.cs, R.files.cs, cs.hh.data.all)
hist(cs.data.all$ESTIMATED.AGE, breaks=25)

rm(cs.hh.data.all)
unique(cs.data.all$COUNTRY)

#unique(paste(cs.data.all$CURRENT.EDU.TYPE.LABEL, cs.data.all$CURRENT.EDU.TYPE))
#unique(paste(cs.data.all$HIGHEST.ATTAINED.EDU.LABEL, cs.data.all$HIGHEST.ATTAINED.EDU))

#temp <- cs.data.all[ cs.data.all$CURRENT.EDU.TYPE == -2 | (as.numeric(cs.data.all$CURRENT.EDU.TYPE) >= as.numeric(cs.data.all$HIGHEST.ATTAINED.EDU) ), ]
#table(cs.data.all$ECON.STATUS.CURR.SELFDEF)

H.files.lon <- c(
  "K:/EU-SILC longitudinal/2 L-2006/UDB_L06H_ver 2006-2 from 01-03-2009.csv",
  "K:/EU-SILC longitudinal/3 L-2007/UDB_l07H_ver 2007-5 from 01-08-2011.csv",
  "K:/EU-SILC longitudinal/4 L-2008/UDB_l08H_ver 2008-4 from 01-03-2012.csv",
  "K:/EU-SILC longitudinal/5 L-2009/UDB_l09H_ver 2009-4 from 01-03-2013.csv",
  "K:/EU-SILC longitudinal/6 L-2010/UDB_l10H_ver 2010-5 from 01-08-2014.csv",
  "K:/EU-SILC longitudinal/7 L-2011/UDB_l11H_ver 2011-4 from 01-03-2015.csv",
  "K:/EU-SILC longitudinal/8 L-2012/UDB_l12H_ver 2012-3 from 01-08-2015.csv",
  "K:/EU-SILC longitudinal/9 L-2013/UDB_l13H_ver 2013-1 from 01-08-2015.csv"
)
lon.hh.data.all <- lifeCourseExposureTrajectories::readHouseholdData(H.files.lon)

unique(lon.hh.data.all$COUNTRY)
lon.hh.data.all <- subset(lon.hh.data.all, COUNTRY %in% c("AT"))
unique(lon.hh.data.all$COUNTRY)

P.files.lon <- c(
  "K:/EU-SILC longitudinal/2 L-2006/UDB_L06P_ver 2006-2 from 01-03-2009.csv",
  "K:/EU-SILC longitudinal/3 L-2007/UDB_l07P_ver 2007-5 from 01-08-2011.csv",
  "K:/EU-SILC longitudinal/4 L-2008/UDB_l08P_ver 2008-4 from 01-03-2012.csv",
  "K:/EU-SILC longitudinal/5 L-2009/UDB_l09P_ver 2009-4 from 01-03-2013.csv",
  "K:/EU-SILC longitudinal/6 L-2010/UDB_l10P_ver 2010-5 from 01-08-2014.csv",
  "K:/EU-SILC longitudinal/7 L-2011/UDB_l11P_ver 2011-4 from 01-03-2015.csv",
  "K:/EU-SILC longitudinal/8 L-2012/UDB_l12P_ver 2012-3 from 01-08-2015.csv",
  "K:/EU-SILC longitudinal/9 L-2013/UDB_l13P_ver 2013-1 from 01-08-2015.csv"
)
data.rshp <- lifeCourseExposureTrajectories::readIndividualSeq(P.files.lon, lon.hh.data.all)
unique(data.rshp$CNTRY)

rm(lon.hh.data.all)

## PL031: Self-defined current economic status
##
## 1 Employee working full-time
## 2 Employee working part-time
## 3 Self-employed working full-time (including family worker)
## 4 Self-employed working part-time (including family worker)
## 5 Unemployed
## 6 Pupil, student, further training, unpaid work experience
## 7 In retirement or in early retirement or has given up business
## 8 Permanently disabled or/and unfit to work
## 9 In compulsory military community or service
## 10 Fulfilling domestic tasks and care responsibilities
## 11 Other inactive person

seqstatl(data.rshp[, 8:11])

# data.alphabet <- seq(1:11)
# data.labels <- c(
#   "Employee working full-time",
#   "Employee working part-time",
#   "Self-employed working full-time (including family worker)",
#   "Self-employed working part-time (including family worker)",
#   "Unemployed",
#   "Pupil, student, further training, unpaid work experience",
#   "In retirement or in early retirement or has given up business",
#   "Permanently disabled or/and unfit to work",
#   "In compulsory military community or service",
#   "Fulfilling domestic tasks and care responsibilities",
#   "Other inactive person"
# )
# data.scodes <- c("EWFT",  "EWPT",  "SEFT",  "SEPT",  "UNEM",  "STUD",  "RETD",  "UNFT",  "CMCS",  "DOME",  "INAC")
data.seq <- seqdef(data.rshp, 8:11, alphabet = data.alphabet[1:11], states = data.scodes[1:11], labels = data.labels[1:11], xtstep = 1)

# markov.chain <- lifeExposureTrajectories::asMarkovChain(data.seq, data.scodes)
# lifeExposureTrajectories::plotMarkovChain(markov.chain)

data.rshp.subset <- lifeCourseExposureTrajectories::stratifiedSeqSample(data.rshp, 0.25)

#help(memory.size)
#memory.limit(size=48000)
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

df <- lifeCourseExposureTrajectories::traversePreOrder(st$root)

edu.data <- lifeCourseExposureTrajectories::getCrossSectEdu(cs.data.all, country.filter = unique(data.rshp.subset$CNTRY))
rm(cs.data.all)
gc()

################################
#
# PARALLEL PART
#
################################



PATH <- "N:\\tfu\\552_HEALS\\Projektarbeit\\WP11"

parallel.lifelong.exposure <- function(i) {
  
  message(individuals$id[i])
  
  cl <- makeCluster(NUM_CORES)
  #cl <- makeCluster( mpi.universe.size() - 2, type="MPI" )
  
  clusterEvalQ(cl, library(lifeCourseExposureTrajectories))
  clusterEvalQ(cl, library(readxl))
  clusterEvalQ(cl, library(parallel))
  
  clusterExport(cl, "NUM_SIM")
  clusterExport(cl, "NUM_CORES")
  clusterExport(cl, "df")
  clusterExport(cl, "st")
  clusterExport(cl, "data.seq")
  clusterExport(cl, "individuals")
  #clusterExport(cl, "par_determineFullTraj")
  #clusterExport(cl, "par_lifeTrajSim")
  
  #sim.results <- NULL
  
  INDIV_SUBJID <- individuals$id[i]
  INDIV_AGE <- individuals$age[i]
  INDIV_SEX <- ifelse(as.character(individuals$sex[i]) == "M", 1, 2)
  
  INDIV_EDU <- individuals$edulevel[i]
  
  econ.stat <- subset(
    edu.data, 
    (ESTIMATED.AGE >= INDIV_AGE-1 & ESTIMATED.AGE <= INDIV_AGE+1) & EDULEVEL == INDIV_EDU & SEX == INDIV_SEX
  )$ECON.STATUS.CURR.SELFDEF
  INDIV_ACT <- data.scodes[ sample(econ.stat, 1) ]
  
  #sim.results.grpd <- par_determineFullTraj(cl, INDIV_SUBJID, INDIV_AGE, INDIV_SEX, INDIV_EDU, INDIV_ACT)
  sim.results.grpd <- lifeCourseExposureTrajectories::determineFullTraj(INDIV_SUBJID, INDIV_AGE, INDIV_SEX, INDIV_EDU, INDIV_ACT)
  
  stopifnot( length(unique(sim.results.grpd$INDIV_SUBJID)) == 1 )
  
  ##write.csv(sim.results.grpd, file = paste0("C:\\TEMP\\sim-out\\lifetraj-", INDIV_SUBJID, ".csv"), row.names=F)
  message("Daily exposures ...")
  
  exposure.all <- lifeCourseExposureTrajectories::getExposureData(
    INDIV_SUBJID, 
    PATH, 
    stressors = c("NO2", "UV", "EMF")
  )
  stopifnot( length(unique(exposure.all$sample_ID)) == 1 )
  
  # impute based on 'surrounding' age for same individual if data for specific age is missing
  
  #sim.lifetraj <- subset(sim.results.grpd, INDIV_SUBJID == INDIV_SUBJID)
  sim.lifetraj <- sim.results.grpd
  
  #clusterExport(cl, "sim.lifetraj")
  #clusterExport(cl, "exposure.all")

  parallel.gapfill.exposure <- function(stressor, INDIV_SUBJID, INDIV_AGE) {
    exposure.add <- data.frame()
    
    for (t in c(1:nrow(sim.lifetraj))) {
      
      sim.age <- sim.lifetraj[ t, ]$AGE
      
      if (sim.age <= INDIV_AGE) {
        sim.emp.scode  <- sim.lifetraj[ t, ]$EMP.scode
        
        exposure.subset <- data.frame()
        age.diff <- 0
        while (nrow(exposure.subset) <= 0 & age.diff < 10) {
          exposure.subset <- subset(
            exposure.all, 
            sample_ID == INDIV_SUBJID & 
              (age >= sim.age - age.diff & age <= sim.age + age.diff) &
              EMP.scode == sim.emp.scode &
              STRESSOR == stressor
          )
          if (age.diff > 0 & nrow(exposure.subset) > 0) {
            exposure.subset$type <- 3
            exposure.subset$comment <- paste0("based on different age (original age=", exposure.subset$age[1], ")")
            exposure.subset$age <- sim.age
            exposure.add <- rbind(exposure.subset, exposure.add) 
          }
          age.diff <- age.diff + 1
        }
      }
      
    }
    return( exposure.add )
  }
  
  exposure.temp <- parLapply(
    cl,
    unique(exposure.all$STRESSOR), 
    parallel.gapfill.exposure, INDIV_SUBJID=INDIV_SUBJID, INDIV_AGE=INDIV_AGE
  )
  t <- do.call(rbind.data.frame, exposure.temp)
  exposure.all <- rbind(exposure.all, t)
  stopifnot( length(unique(exposure.all$sample_ID)) == 1 )
  
  #sort(unique(exposure.all[ exposure.all$sample_ID == INDIV_SUBJID & exposure.all$STRESSOR == sim.stressor, ]$age))
  #summary(exposure.all)
  
  #
  sim.exposure.all <- merge(
    x = sim.results.grpd,
    y = exposure.all,
    by.x = c("INDIV_SUBJID", "AGE", "EMP.scode"),
    by.y = c("sample_ID", "age", "EMP.scode")
    #suffixes = c("SIM", "EXP"),
  )
  stopifnot( length(unique(sim.exposure.all$INDIV_SUBJID)) == 1 )
  #sort(unique(sim.exposure.all$AGE))
  
  library(plyr)
  
  # determine relative weight of exposure estimate per (age, EMP.scode)-combination based on
  # the number of original diaries that were used for estimation
  #
  # equal weight could be reached by replacing the 'aggregate' statement below by
  #   y = count(sim.exposure.PM25, vars = c("INDIV_SUBJID", "AGE", "EMP.scode"))
  #
  sim.exposure.all <- merge(
    x = sim.exposure.all,
    y = aggregate(count ~ INDIV_SUBJID + STRESSOR + AGE + EMP.scode, data=sim.exposure.all[,c("INDIV_SUBJID", "STRESSOR", "AGE", "EMP.scode", "count")], sum),
    by = c("INDIV_SUBJID", "STRESSOR", "AGE", "EMP.scode"),
    suffixes = c(".per.emp", ".all")
  )
  sim.exposure.all$EMP.map.prob <- sim.exposure.all$count.per.emp / sim.exposure.all$count.all
  
  # determine total probability per row ...
  sim.exposure.all$total.prob <- sim.exposure.all$EMP.probability * sim.exposure.all$EMP.map.prob
  
  ##write.csv(sim.exposure.all, file = paste0("C:\\TEMP\\sim-out\\exposure-map-", INDIV_SUBJID, ".csv"), row.names=F)
  
  # ... and sample accordingly
  d <- unique(sim.exposure.all[ c("INDIV_SUBJID", "STRESSOR", "AGE") ])
  d <- d[ with(d, order(STRESSOR, -AGE)), ]
  rownames(d) <- NULL
  
  d$EXP_2.5PCT <- -1
  d$EXP_25PCT <- -1
  d$EXP_MEDIAN <- -1
  d$EXP_75PCT <- -1
  d$EXP_97.5PCT <- -1
  d$EXP_MEAN <- -1
  d$EXP_SD <- -1
  
  library(triangle)
  
  sample.exp.all <- data.frame(INDIV_SUBJID = NULL, STRESSOR = NULL, value = NULL, age = NULL)
  
  parallel.sample.exposure <- function(j) {
    subjid   <- d[ j, ]$INDIV_SUBJID
    stressor <- d[ j, ]$STRESSOR
    age      <- d[ j, ]$AGE
    
    sim.exposure.subset <- subset(sim.exposure.all, INDIV_SUBJID == subjid & STRESSOR == stressor & AGE == age)
    
    sample.rownames <- sample(
      #size = 200,
      size = 100,
      x = rownames(sim.exposure.subset),
      prob = sim.exposure.subset$total.prob,
      replace = T
    )
    
    sample.exp <- data.frame(value = NULL)
    parallel.sample.exposure.inner <- function(r) {
      return(
        data.frame(
          rtriangle(
            #n = 400,
            n = 150,
            a = min(sim.exposure.subset[ r, ]$X2.5._percentile, sim.exposure.subset[ r, ]$mean),
            b = max(sim.exposure.subset[ r, ]$mean, sim.exposure.subset[ r, ]$X97.5._percentile),
            c = median(c(sim.exposure.subset[ r, ]$X2.5._percentile, sim.exposure.subset[ r, ]$mean, sim.exposure.subset[ r, ]$X97.5._percentile))
          )
        )
      )
    }
    sample.exp.temp <- mclapply(
      sample.rownames,
      FUN=parallel.sample.exposure.inner
    )
    sample.exp <- do.call(rbind.data.frame, sample.exp.temp)
    names(sample.exp) <- c("value")
    
    d[ j, c("EXP_2.5PCT", "EXP_25PCT", "EXP_MEDIAN", "EXP_75PCT", "EXP_97.5PCT", "EXP_MEAN", "EXP_SD") ] <- c(quantile(sample.exp$value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(sample.exp$value), sd(sample.exp$value))
    
    return(
      data.frame(
        INDIV_SUBJID = subjid, 
        STRESSOR = stressor, 
        value = sample.exp, 
        age = age
      )
    )
  }
  
  clusterEvalQ(cl, library(triangle))
  clusterEvalQ(cl, library(parallel))
  #clusterExport(cl, "d")
  #clusterExport(cl, "sim.exposure.all")
  
  sample.exp.t <- parLapply(
    cl,
    c(1:nrow(d)),
    parallel.sample.exposure
  )
  t <- do.call(rbind.data.frame, sample.exp.t)
  sample.exp.all <- rbind(sample.exp.all, t)
  
  #names(d)
  #sort(unique(sample.exp.all$age))
  #summary(sample.exp.all)
  
  ##write.csv(sample.exp.all, file = paste0("N:\\tfu\\552_HEALS\\Projektarbeit\\WP11\\LET\\sim-out\\exposure-samples-", INDIV_SUBJID, ".csv"), row.names=F)
  
  # name critical life stages ...
  sample.exp.all$CRITICAL_LIFE_STAGE <- ""
  
  if (nrow(sample.exp.all[ sample.exp.all$age <= 3, ]) > 0) {
    sample.exp.all[ sample.exp.all$age <= 3, ]$CRITICAL_LIFE_STAGE <- "Infancy (1-3)"
  }
  
  if (nrow(sample.exp.all[ sample.exp.all$age > 3 & sample.exp.all$age <= 11, ]) > 0) {
    sample.exp.all[ sample.exp.all$age > 3 & sample.exp.all$age <= 11, ]$CRITICAL_LIFE_STAGE <- "Childhood (4-11)"
  }
  
  if (nrow(sample.exp.all[ sample.exp.all$age > 11 & sample.exp.all$age <= 17, ]) > 0) {
    sample.exp.all[ sample.exp.all$age > 11 & sample.exp.all$age <= 17, ]$CRITICAL_LIFE_STAGE <- "Adolescence (12-17)"
  }
  
  if (nrow(sample.exp.all[ sample.exp.all$age > 17 & sample.exp.all$age <= 39, ]) > 0) {
    sample.exp.all[ sample.exp.all$age > 17 & sample.exp.all$age <= 39, ]$CRITICAL_LIFE_STAGE <- "Adulthood before 40 (18-39)"
  }
  
  if (nrow(sample.exp.all[ sample.exp.all$age > 39 & sample.exp.all$age <= 64, ]) > 0) {
    sample.exp.all[ sample.exp.all$age > 39 & sample.exp.all$age <= 64, ]$CRITICAL_LIFE_STAGE <- "Adulthood before 65 (40-64)"
  }
  
  if (nrow(sample.exp.all[ sample.exp.all$age > 64, ]) > 0) {
    sample.exp.all[ sample.exp.all$age > 64, ]$CRITICAL_LIFE_STAGE <- "Adulthood 65 and older (>=65)"
  }
  
  stopCluster(cl)
  
  # ... and generate aggregate stats
  sample.exp.stats <- aggregate(value ~ INDIV_SUBJID + STRESSOR + CRITICAL_LIFE_STAGE, data=sample.exp.all, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  sample.exp.stats <- cbind(sample.exp.stats[, 1:3 ], data.frame(sample.exp.stats$value))
  
  sample.exp.mean <- aggregate(value ~ INDIV_SUBJID + STRESSOR + CRITICAL_LIFE_STAGE, data=sample.exp.all, FUN=mean)
  names(sample.exp.mean$value) <- "value.MEAN"
  
  sample.exp.sd <- aggregate(value ~ INDIV_SUBJID + STRESSOR + CRITICAL_LIFE_STAGE, data=sample.exp.all, FUN=sd)
  names(sample.exp.sd$value) <- "value.SD"
  
  sample.exp.stats <- merge(
    sample.exp.stats,
    sample.exp.mean,
    by = c("INDIV_SUBJID", "STRESSOR", "CRITICAL_LIFE_STAGE")
  )
  sample.exp.stats <- merge(
    sample.exp.stats,
    sample.exp.sd,
    by = c("INDIV_SUBJID", "STRESSOR", "CRITICAL_LIFE_STAGE")
  )
  
  names(sample.exp.stats) <- c("INDIV_SUBJID", "STRESSOR", "CRITICAL_LIFE_STAGE", "value.2.5PCT", "value.25PCT", "value.50PCT", "value.75PCT", "value.97.5PCT", "value.MEAN", "value.SD")
  ##write.csv(sample.exp.stats, file = paste0("C:\\TEMP\\sim-out\\exposure-stats-", INDIV_SUBJID, ".csv"), row.names=F)
  
  #return( sample.exp.all )
  return( sample.exp.stats )
}

individuals <- read.csv(paste0(PATH, "\\Stream 5 data\\stream5SampleData3.csv"))

library(readxl)

NUM_SIM <- 100

library(parallel)

NUM_MASTER_CORES <- 2
NUM_CORES <- parallel::detectCores() - NUM_MASTER_CORES - 1

cl.master <- makeCluster(NUM_MASTER_CORES)

clusterEvalQ(cl.master, library(lifeCourseExposureTrajectories))
clusterEvalQ(cl.master, library(readxl))
clusterEvalQ(cl.master, library(parallel))

clusterExport(cl.master, "NUM_SIM")
clusterExport(cl.master, "NUM_CORES")
clusterExport(cl.master, "PATH")
clusterExport(cl.master, "df")
clusterExport(cl.master, "st")
clusterExport(cl.master, "data.seq")
clusterExport(cl.master, "edu.data")
clusterExport(cl.master, "individuals")

require(tictoc)
tic()

sim.results <- parLapply(
  cl.master,
  402:403, 
  parallel.lifelong.exposure
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

