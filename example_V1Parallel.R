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
#
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
if(!all(c('MGDrivE','MGDrivE2') %in% installed.packages()[ ,'Package']) ){
  stop("Packages 'MGDrivE', 'MGDrivE2', and are required.")
}


###############################################################################
## Experimental Setup and Paths Definition
###############################################################################
USER <- 1
numRep <- 5
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
simTime <- 1*365
sampTime <- 1

# gd releases
relStart <- 60
relInt <- 7

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
                                    'mmrF' = 0.05,
                                    'numRel' = seq.int(from = 0, to = 5, by = 1),
                                    'sizeRel' = c(0.1)* totPopSize ))

numPC <- NROW(paramCombo)

# 46min


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
                                     formatC(x = paramCombo[[x,'numRel']], width = 2, format = 'd', flag = '0'),
                                     formatC(x = paramCombo[[x,'sizeRel']], width = 4, format = 'd', flag = '0')),
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
    if(paramCombo[x,'numRel'] !=0){
      patchReleases[[1]]$maleReleases <- MGDrivE::generateReleaseVector(driveCube = cube,
                                                                        releasesParameters = list('releasesStart'=relStart,
                                                                                                  'releasesNumber'=paramCombo[x,'numRel'],
                                                                                                  'releasesInterval'=relInt,
                                                                                                  'releaseProportion'=paramCombo[x,'sizeRel']))
    }

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
# what kinds of analysis here?
#  function to reduce by carrier/non-carrier?
#  functions to reduce by allele-count?
#  separate or combine by sex
#  Remove 0 alleles?
#    Is this covered by the reductions? Do the reductions need to count all alleles?
#    If they don't, what happens to the statistics?







########################################
### Plotting
########################################
# something that works for now
# redo with ggplot2
# Throw save option in here - steal from fitting scripts?

# get dirs
allDirs <- list.dirs(path = outDir, full.names = TRUE, recursive = FALSE)

# plot all reps from first parameter set
# MGDrivE::plotMGDrivEMult(readDir = allDirs[1], lwd = 0.35, alpha = 0.75)











