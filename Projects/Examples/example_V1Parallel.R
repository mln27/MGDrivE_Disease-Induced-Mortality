###############################################################################
#   _    _____      ______                           __
#  | |  / <  /     / ____/  ______ _____ ___  ____  / /__
#  | | / // /_____/ __/ | |/_/ __ `/ __ `__ \/ __ \/ / _ \
#  | |/ // /_____/ /____>  </ /_/ / / / / / / /_/ / /  __/
#  |___//_/     /_____/_/|_|\__,_/_/ /_/ /_/ .___/_/\___/
#                                         /_/
# .............................................................................
# V1 Example Script
# Ndeffo Mbah Lab
# Jared Bennett, jared_bennett@berkeley.edu
#
# 20220412
# Initial setup.
# This script is a parallel implementation of V1 for testing the SEM design.
#
# 20220424
# Switched the cube - V1 will work for all of the trans constructs, but not the cis.
# Checked releases, added basic plotting at the bottom for initial testing.
# The non-reduced cube takes several minutes to run.
#
# 20220526
# I've made a bunch of changes, and apparently didn't take notes.
# Latest addition is allele counts in the analysis.
#
#
###############################################################################
## Setup packages
###############################################################################
rm(list=ls()); gc()
startTime <- Sys.time()
set.seed(1029384756)

# Check that required packages are installed
# This does not check that it's the proper version of V2!!!!
if(!all(c('MGDrivE','MGDrivE2','ggplot2') %in% installed.packages()[ ,'Package']) ){
  stop("Packages 'MGDrivE', 'MGDrivE2', and 'ggplot2' are required.")
}


###############################################################################
## Experimental Setup and Paths Definition
###############################################################################
USER <- 1
numRep <- 3
numCores <- 2

###############################################################################
if(USER == 1){
  # Jared Laptop
  baseOut<-'~/Desktop/OUTPUT/MGDrivE/V1_3'
}else if(USER == 2){
  # ??
  baseOut<-''
}else{
  warning('invalid user!')
}

repNames <- formatC(x = 1:numRep, width = 3, format = 'd', flag = '0')


###############################################################################
## Inputs
###############################################################################
########################################
## Raw Inputs
########################################
simTime <- 2*365
sampTime <- 1

# gd releases
relStart1 <- 60
relInt1 <- 7

# sem releases
relStart2 <- 300
relInt2 <- 7

# default bio parameters and population size
bioParameters <- list('betaK'=20, 'tEgg'=5, 'tLarva'=6, 'tPupa'=4, 'popGrowth'=1.175, 'muAd'=0.09)
totPopSize <- 10000

# sweep over migration rates - these are adult lifetime rates
migRates <- c(0.01)
numPatch <- 3
patchSet <- 1:(numPatch-1)

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- MGDrivE::basicBatchMigration(batchProbs=0,
                                               sexProbs=c(.5,.5),
                                               numPatches=numPatch)


########################################
## Sweeps
########################################
# Setup desired parameter ranges here
#  This then creates a dataframe with every combination

# sweep over cube changing parameters
#  pF:
#  qF:
#  rF:
#  aF:
#  bF:
#  cF:
#  xF:
#  yF:
#  mmrF:
#  numRel: number of releases
#  sizeRel: size of releases, as a percentage of total population * totPopSize
paramCombo <- as.matrix(expand.grid('pF' = c(0.9),
                                    'qF' = c(0.9),
                                    'rF' = c(0.05),
                                    'aF' = c(0.9),
                                    'bF' = c(0.9),
                                    'cF' = c(0.9),
                                    'xF' = c(0.9),
                                    'yF' = c(0.9),
                                    'mmrF' = 0.00,
                                    'numRel1' = seq.int(from = 0, to = 2, by = 1),
                                    'sizeRel1' = c(0.1)* totPopSize,
                                    'numRel2' = seq.int(from = 0, to = 2, by = 1),
                                    'sizeRel2' = c(0.1)* totPopSize ))

numPC <- NROW(paramCombo)


###############################################################################
### Simulation
###############################################################################
########################################
### Migration
########################################
for(mR in migRates){
  # get initial time
  startTimeLand <- Sys.time()

  # moveMat
  if(numPatch==1){
    # single-patch
    moveMat <- matrix(data = 1, nrow = numPatch, ncol = numPatch)
  } else {
    # multiple patches, connected in a line
    # convert lifetime rate to daily rate
    dmR <- 1-(1-mR)^(bioParameters$muAd)
    # setup movement matrix
    moveMat <- matrix(data = 0, nrow = numPatch, ncol = numPatch)
    moveMat[cbind(patchSet, patchSet+1)] <- dmR
    moveMat[cbind(patchSet+1, patchSet)] <- dmR
    diag(moveMat) <- 1 - rowSums(moveMat)
  }

  # build output directories names
  mRName <- formatC(x = mR*totPopSize, width = 3, format = 'd', flag = '0')
  outDir <- file.path(baseOut, mRName)

  # create base directory
  if(!dir.exists(outDir)){ dir.create(path = outDir, recursive = TRUE) }


  ########################################
  ### Run Experiments
  ########################################
  # Setup the cluster
  #  have to restart each time to validate changed objects
  #  https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
  if(Sys.info()[['sysname']] == "Windows"){
    # windows can't run fork clusters
    #  it has to use sockets
    cl <- parallel::makePSOCKcluster(names = numCores)

    # load libraries on each socket
    #  if we call using package::function(), this is unnecessary
    # parallel::clusterEvalQ(cl = cl, expr = {library(MGDrivE2)})

    # export required objects to each socket
    parallel::clusterExport(
      cl = cl,
      varlist = c("outDir", "paramCombo","repNames",
                  "moveMat", "relStart", "relInt",
                  "simTime", "sampTime", "bioParameters", "totPopSize",
                  "batchMigration", 'numRep')
    )

  } else {
    # *nix systems
    # Sys.info() does distuinguish Linux vs Mac, but we don't care for clusters
    # no copying, no memory issues!
    cl <- parallel::makeForkCluster(nnodes = numCores)

  }

  # Set parallel seed
  parallel::clusterSetRNGStream(cl = cl, iseed = sample(x = (-.Machine$integer.max):.Machine$integer.max,
                                                        size = 1, replace = FALSE))

  # Run
  retList <- parallel::clusterApplyLB(cl = cl, x = 1:numPC, fun = function(x){

    ####################
    # Folders
    ####################
    # generate experiment directory with all of the simulation parameters
    # This ensures that the directories are unique
    # made sure the width is enough that values are not cut off
    runFolders <- file.path(outDir,
                            paste0(c(formatC(x = paramCombo[x,c('pF','qF','rF','aF','bF','cF','xF','yF','mmrF')]*1000,
                                             width = 4, format = 'd', flag = '0'),
                                     formatC(x = paramCombo[[x,'numRel1']], width = 2, format = 'd', flag = '0'),
                                     formatC(x = paramCombo[[x,'sizeRel1']], width = 4, format = 'd', flag = '0'),
                                     formatC(x = paramCombo[[x,'numRel2']], width = 2, format = 'd', flag = '0'),
                                     formatC(x = paramCombo[[x,'sizeRel2']], width = 4, format = 'd', flag = '0')),
                                   collapse = '_'),
                            repNames)

    # make folders
    for(i in runFolders){ dir.create(path = i, recursive = TRUE) }

    ####################
    # Build cube
    ####################
    # equal parameters between males and females
    #  all deposition is turned off
    #  no genotype-specific fitness costs
    cube <- MGDrivE2::cubeSEMtrans(pF = paramCombo[[x,'pF']], qF = paramCombo[[x,'qF']], rF = paramCombo[[x,'rF']],
                                   aF = paramCombo[[x,'aF']], bF = paramCombo[[x,'bF']], cF = paramCombo[[x,'cF']],
                                   xF = paramCombo[[x,'xF']], yF = paramCombo[[x,'yF']],
                                   mmrF = paramCombo[[x,'mmrF']])

    ####################
    # Releases
    ####################
    # default release list required - set all values to NULL
    patchReleases <- replicate(n=NROW(moveMat),
                               expr={list(maleReleases=NULL, femaleReleases=NULL,
                                          eggReleases=NULL, matedFemaleReleases=NULL,
                                          iSpray=NULL, smSpray=NULL)},
                               simplify=FALSE )

    # generate releases at the specified times
    # safety for when there are no releases
    if(paramCombo[x,'numRel1'] !=0){
      patchReleases[[1]]$maleReleases <- MGDrivE::generateReleaseVector(driveCube = cube,
                                                                        releasesParameters = list('releasesStart'=relStart1,
                                                                                                  'releasesNumber'=paramCombo[x,'numRel1'],
                                                                                                  'releasesInterval'=relInt1,
                                                                                                  'releaseProportion'=paramCombo[x,'sizeRel1']))

      # check for SEM releases
      #  if we don't do GD releases, no point in doing SEM releases
      if(paramCombo[x,'numRel2'] !=0){
        # again, do male releases

        # swap the "releaseType"in the cube for now, so that we don't have to mess with
        #  the underlying code
        oRT <- cube$releaseType
        cube$releaseType <- "WWHH"

        patchReleases[[1]]$maleReleases <- c(patchReleases[[1]]$maleReleases,
                                             MGDrivE::generateReleaseVector(driveCube = cube,
                                                                            releasesParameters = list('releasesStart'=relStart2,
                                                                                                      'releasesNumber'=paramCombo[x,'numRel2'],
                                                                                                      'releasesInterval'=relInt2,
                                                                                                      'releaseProportion'=paramCombo[x,'sizeRel2'])))
        # reset release type
        cube$releaseType <- oRT

      } # end SEM releases
    } # end GD releases

    ####################
    # Parameters
    ####################
    netPar <- MGDrivE::parameterizeMGDrivE(runID=1, simTime=simTime, sampTime = sampTime, nPatch=NROW(moveMat),
                                           beta=bioParameters$betaK, muAd=bioParameters$muAd,
                                           popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                                           tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                                           AdPopEQ=totPopSize, inheritanceCube = cube)

    ####################
    # Run
    ####################
    # set MGDrivE to run stochastic
    MGDrivE::setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

    # build network
    MGDrivESim <- MGDrivE::Network$new(params=netPar,
                                       driveCube=cube,
                                       patchReleases=patchReleases,
                                       migrationMale=moveMat,
                                       migrationFemale=moveMat,
                                       migrationBatch=batchMigration,
                                       directory=runFolders,
                                       verbose = FALSE)
    # run
    MGDrivESim$multRun(verbose = FALSE)

    ####################
    # Analyze
    ####################
    # First, split output by patch
    # Second, aggregate females by their mate choice
    for(i in 1:numRep){
      MGDrivE::splitOutput(readDir = runFolders[i], remFile = TRUE, verbose = FALSE)
      MGDrivE::aggregateFemales(readDir = runFolders[i], genotypes = cube$genotypesID,
                                remFile = TRUE, verbose = FALSE)
    }

    gc()

  }) # end cluster apply

  # stop cluster
  parallel::stopCluster(cl)

  # print time it took
  print(paste0('  Sims: ', capture.output(difftime(time1 = Sys.time(), time2 = startTimeLand))))

}# end loop over migration rates

# print time it took
print(paste0('Total: ', capture.output(difftime(time1 = Sys.time(), time2 = startTime))))


###############################################################################
### Analysis and Plotting
###############################################################################
########################################
### Analysis
########################################
####################
# Data
####################
# grab the parameter directories
migDirs <- list.dirs(path = baseOut, full.names = TRUE, recursive = FALSE)
migParamDirs <- lapply(X = migDirs, FUN = list.dirs, full.names = TRUE, recursive = FALSE)

####################
# Genetic Stuff
####################
# taken from the cube, so I didn't have to rerun it
genos <- c(
"WWWW","GWWW","UWWW","RWWW","VWWW","SWWW","GGWW","GUWW","GRWW","GVWW",
"GSWW","UUWW","RUWW","UVWW","SUWW","RRWW","RVWW","RSWW","VVWW","SVWW",
"SSWW","WWHW","GWHW","UWHW","RWHW","VWHW","SWHW","GGHW","GUHW","GRHW",
"GVHW","GSHW","UUHW","RUHW","UVHW","SUHW","RRHW","RVHW","RSHW","VVHW",
"SVHW","SSHW","WWRW","GWRW","UWRW","RWRW","VWRW","SWRW","GGRW","GURW",
"GRRW","GVRW","GSRW","UURW","RURW","UVRW","SURW","RRRW","RVRW","RSRW",
"VVRW","SVRW","SSRW","WWEW","GWEW","UWEW","RWEW","VWEW","SWEW","GGEW",
"GUEW","GREW","GVEW","GSEW","UUEW","RUEW","UVEW","SUEW","RREW","RVEW",
"RSEW","VVEW","SVEW","SSEW","WWHH","GWHH","UWHH","RWHH","VWHH","SWHH",
"GGHH","GUHH","GRHH","GVHH","GSHH","UUHH","RUHH","UVHH","SUHH","RRHH",
"RVHH","RSHH","VVHH","SVHH","SSHH","WWHR","GWHR","UWHR","RWHR","VWHR",
"SWHR","GGHR","GUHR","GRHR","GVHR","GSHR","UUHR","RUHR","UVHR","SUHR",
"RRHR","RVHR","RSHR","VVHR","SVHR","SSHR","WWEH","GWEH","UWEH","RWEH",
"VWEH","SWEH","GGEH","GUEH","GREH","GVEH","GSEH","UUEH","RUEH","UVEH",
"SUEH","RREH","RVEH","RSEH","VVEH","SVEH","SSEH","WWRR","GWRR","UWRR",
"RWRR","VWRR","SWRR","GGRR","GURR","GRRR","GVRR","GSRR","UURR","RURR",
"UVRR","SURR","RRRR","RVRR","RSRR","VVRR","SVRR","SSRR","WWER","GWER",
"UWER","RWER","VWER","SWER","GGER","GUER","GRER","GVER","GSER","UUER",
"RUER","UVER","SUER","RRER","RVER","RSER","VVER","SVER","SSER","WWEE",
"GWEE","UWEE","RWEE","VWEE","SWEE","GGEE","GUEE","GREE","GVEE","GSEE",
"UUEE","RUEE","UVEE","SUEE","RREE","RVEE","RSEE","VVEE","SVEE","SSEE")

# build the genotypes of interest
#  have to do this because there's no total by default, and I want to see that
goiList <- c(list("Total"=genos),setNames(object = as.list(genos), nm = genos)
             )

# build alleles of interest
#  These functions pull out all combinations of alleles.
#  we'll use them to create the goiLists for other analysis
#  Could probably put these in a loop and generate the list that way. Oh well.
one_W <- c(grep(pattern = "W...|.W..", x = genos, value = TRUE),
           grep(pattern = "WW..", x = genos, value = TRUE))

one_G <- c(grep(pattern = "G...|.G..", x = genos, value = TRUE),
           grep(pattern = "GG..", x = genos, value = TRUE))

one_U <- c(grep(pattern = "U...|.U..", x = genos, value = TRUE),
           grep(pattern = "UU..", x = genos, value = TRUE))

one_R <- c(grep(pattern = "R...|.R..", x = genos, value = TRUE),
           grep(pattern = "RR..", x = genos, value = TRUE))

one_V <- c(grep(pattern = "V...|.V..", x = genos, value = TRUE),
           grep(pattern = "VV..", x = genos, value = TRUE))

one_S <- c(grep(pattern = "S...|.S..", x = genos, value = TRUE),
           grep(pattern = "SS..", x = genos, value = TRUE))


two_W <- c(grep(pattern = "..W.|...W", x = genos, value = TRUE),
           grep(pattern = "..WW", x = genos, value = TRUE))

two_H <- c(grep(pattern = "..H.|...H", x = genos, value = TRUE),
           grep(pattern = "..HH", x = genos, value = TRUE))

two_R <- c(grep(pattern = "..R.|...R", x = genos, value = TRUE),
           grep(pattern = "..RR", x = genos, value = TRUE))

two_E <- c(grep(pattern = "..E.|...E", x = genos, value = TRUE),
           grep(pattern = "..EE", x = genos, value = TRUE))


goiList_1 <- list("W" = one_W,
                  "G" = one_G,
                  "U" = one_U,
                  "R" = one_R,
                  "V" = one_V,
                  "S" = one_S,
                  "Total" = c(one_W, one_G, one_U, one_R, one_V, one_S))

goiList_2 <- list("W" = two_W,
                  "H" = two_H,
                  "R" = two_R,
                  "E" = two_E,
                  "Total" = c(two_W, two_H, two_R, two_E))

goiList_both <- list("one_W" = one_W,
                     "one_G" = one_G,
                     "one_U" = one_U,
                     "one_R" = one_R,
                     "one_V" = one_V,
                     "one_S" = one_S,
                     "two_W" = two_W,
                     "two_H" = two_H,
                     "two_R" = two_R,
                     "two_E" = two_E,
                     "Total" = c(two_W, two_H, two_R, two_E))
# the "Total" here is kinda a cheat
#  since we always have a locus 1 paired with a locus 2, the total counts for both
#  have to be equal. I used the total from two cause it was less typing.

####################
# Analysis
####################


# loop over parameter directories, do analysis!
#  coudl also loop this better - I'm lazy.
for(wDir in migParamDirs[[1]]){
  # genotype analysis
  MGDrivE2::analyze_ggplot_CSV(read_dir = wDir, name = "genos", sex = "both",
                               patch_agg = FALSE, goi = goiList, drop_zero_goi = TRUE)

  # locus 1
  MGDrivE2::analyze_ggplot_CSV(read_dir = wDir, name = "allele_1", sex = "both",
                               patch_agg = FALSE, goi = goiList_1, drop_zero_goi = TRUE)

  # locus 2
  MGDrivE2::analyze_ggplot_CSV(read_dir = wDir, name = "allele_2", sex = "both",
                               patch_agg = FALSE, goi = goiList_2, drop_zero_goi = TRUE)

  # both loci
  MGDrivE2::analyze_ggplot_CSV(read_dir = wDir, name = "allele_both", sex = "both",
                               patch_agg = FALSE, goi = goiList_both, drop_zero_goi = TRUE)
}







########################################
### Plotting
########################################
library(ggplot2)

####################
# Data
####################
# get dirs
migDirs <- list.dirs(path = baseOut, full.names = TRUE, recursive = FALSE)
migAnalysisDirs <- lapply(X = migDirs, FUN = function(x){
  grep(pattern = "analysis_", x = list.dirs(path = x, full.names = TRUE), value = TRUE, fixed = TRUE)
})

# read 2 sets of summary data
# columns: "Repetitions" "Sex" "Patch" "Time" "GOI" "Count"
gPDat <- read.csv(file = file.path(migAnalysisDirs[[1]][3], "plot.csv"),
                  header = TRUE, sep = ",")
# columns: "Sex" "Patch" "Time" "GOI" "Min" "Mean" "Max" "2.5%" "50%" "97.5%"
gMDat <- read.csv(file = file.path(migAnalysisDirs[[1]][3], "metrics.csv"),
                  header = TRUE, sep = ",", check.names = FALSE)

####################
# Traces
####################
# data
ggplot(data = gPDat ) +
# color scheme and organization
geom_path(aes(x = Time, y = Count, group = interaction(Repetitions, GOI),
              color = GOI), alpha = 0.35) +
facet_grid(Patch ~ Sex, scales = "free_y") +

# # mean
# geom_line(data = gMDat, mapping = aes(x = Time, y = Mean,
#                                       group = interaction(GOI),
#                                       color = GOI)) +

# aesthetics
theme_bw(base_size = 16) +
ggtitle("Population Dynamics")


####################
# Mean + Quantiles
####################
# https://stackoverflow.com/questions/28586635/shade-region-between-two-lines-with-ggplot
# https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot

# data + color scheme and organization
ggplot(data = gMDat, aes(x = Time, y = Mean, group = interaction(GOI),
              color = GOI), alpha = 0.15 ) +
facet_grid(Patch ~ Sex, scales = "free_y") +

# min/max lines
#  not necessary, as geom_ribbon plots them as well
# geom_line(data = gMDat, mapping = aes(x = Time, y = get("2.5%"),
#                                       group = interaction(GOI),
#                                       color = GOI), size = 0.5) +
# geom_line(data = gMDat, mapping = aes(x = Time, y = get("97.5%"),
#                                       group = interaction(GOI),
#                                       color = GOI)) +

# ribbon and fill
geom_ribbon(data = gMDat, aes(ymin = get("2.5%"), ymax = get("97.5%"), fill = GOI),
            alpha = 0.25, size = 0.0) +

# mean
geom_line(data = gMDat, mapping = aes(x = Time, y = Mean,
                                      group = interaction(GOI),
                                      color = GOI), size = 0.75) +
# aesthetics
theme_bw(base_size = 16) +
ggtitle("Population Dynamics")



