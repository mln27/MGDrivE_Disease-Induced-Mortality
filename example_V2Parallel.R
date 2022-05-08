###############################################################################
#   _    _____        ______                           __
#  | |  / /__ \      / ____/  ______ _____ ___  ____  / /__
#  | | / /__/ /_____/ __/ | |/_/ __ `/ __ `__ \/ __ \/ / _ \
#  | |/ // __/_____/ /____>  </ /_/ / / / / / / /_/ / /  __/
#  |___//____/    /_____/_/|_|\__,_/_/ /_/ /_/ .___/_/\___/
#                                           /_/
# .............................................................................
# V2 Example Script
# Ndeffo Mbah Lab
# Jared Bennett, jared_bennett@berkeley.edu
#
# 20220328
# Initial setup.
#
# 20220424
# Adapt for parallel execution. This was stolen from the SEM test script and then
# modeled after the V1 script.
# However, it runs parallel over reps, not parameter sets.
#
# 20220425
# Finished the adaptation. Setup a more flexible method for doing SEM "releases".
# This will still require significant manual intervention.
# Works on windows and unix, need to get analysis and plotting stuff going.
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
if(!all(c('MGDrivE','MGDrivE2','ggplot2','Matrix') %in% installed.packages()[ ,'Package']) ){
  stop("Packages 'MGDrivE', 'MGDrivE2', 'Matrix', and 'ggplot2' are required.")
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
  baseOut<-'~/Desktop/OUTPUT/MGDrivE/V2'
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
dtStoch <- 0.1

# gd releases
relStart1 <- 60
relInt1 <- 7

# sem spraying
relStart2 <- 300
relInt2 <- 1

# default bio parameters and population size
bioParameters <- list('betaK'=20, 'tEgg'=5, 'tLarva'=6, 'tPupa'=4, 'popGrowth'=1.175, 'muAd'=0.09)
numFem <- 5000 # adult females is 1/2 the total population size

# shape parameters
# variation in stages - defines how many compartments in the Erlang distribution
#  Some tweaking of these parameters may help convergence with V1
cvStage <- c('E'=0.6, 'L'=0.4, 'P'=0.7)
shapeStage <- MGDrivE2::get_shape(cv = cvStage,
                                  q = 1/unlist(bioParameters[c('tEgg','tLarva','tPupa')], use.names = FALSE))

# build parameter list
bioPars <- list(
  qE = 1/bioParameters$tEgg,
  nE = shapeStage[1],
  qL = 1/bioParameters$tLarva,
  nL = shapeStage[2],
  qP = 1/bioParameters$tPupa,
  nP = shapeStage[3],
  beta = bioParameters$betaK,
  muF = bioParameters$muAd,
  muM = bioParameters$muAd,
  nu = 1/(4/24),
  phi = 0.5
)

# calc density-independent aquatic death, based on adult parameters
bioPars$muE <- bioPars$muL <- bioPars$muP <- MGDrivE2::solve_muAqua(params = bioPars,
                                                                    rm = bioParameters$popGrowth)

# sweep over migration rates - these are adult lifetime rates
migRates <- c(0.01)
numPatch <- 3
patchSet <- 1:(numPatch-1)

# this is specific to the inheritance pattern!!!
#  the first row is the "off" genotypes - just GD, no SEM
#  the second row is what those genotypes become when turned "on" - GE and SEM
#  these have to be in the genotypes of the cube being simulated
swapMat <- matrix(data = c("GW","GG","GU","GR","GV","GH","GS",
                           "HW","HH","HU","HR","HV","HH","HS"),
                  nrow = 7, ncol = 2, dimnames = list(NULL, c("off","on")))


########################################
## Sweeps
########################################
# Setup desired parameter ranges here
#  This then creates a dataframe with every combination
#  This may be dependent on the inheritance pattern as well

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
                                    'mmrF' = 0.05,
                                    'numRel1' = seq.int(from = 0, to = 5, by = 1),
                                    'sizeRel1' = c(0.1)* numFem*2,
                                    'numRel2' = seq.int(from = 30, to = 60, by = 30) ))

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

  ####################
  # Movement Setup
  ####################
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

  # convert probability matrix to rate matrix
  probRateList <- MGDrivE2::movement_prob2rate(tau = moveMat)
  # grab adjacency matrix to setup transitions
  adjMat <- as(probRateList$mat, 'ngCMatrix')


  # build output directories names
  mRName <- formatC(x = mR*numFem*2, width = 3, format = 'd', flag = '0')
  outDir <- file.path(baseOut, mRName)

  # create base directory
  if(!dir.exists(outDir)){ dir.create(path = outDir, recursive = TRUE) }


  ########################################
  ### Parameter Sweeps
  ########################################
  for(nP in 1:numPC){
    # get initial time
    startTimeParams <- Sys.time()

    ####################
    # Folders
    ####################
    # generate experiment directory with all of the simulation parameters
    # This ensures that the directories are unique
    # made sure the width is enough that values are not cut off
    runDir <- file.path(outDir,
                        paste0(c(formatC(x = paramCombo[nP,c('pF','qF','rF','aF','bF','cF','mmrF')]*1000,
                                         width = 4, format = 'd', flag = '0'),
                                 formatC(x = paramCombo[[nP,'numRel1']], width = 2, format = 'd', flag = '0'),
                                 formatC(x = paramCombo[[nP,'sizeRel1']], width = 4, format = 'd', flag = '0'),
                                 formatC(x = paramCombo[[nP,'numRel2']], width = 4, format = 'd', flag = '0')),
                               collapse = '_'))

    runFolders <- file.path(runDir, repNames)

    # make folders
    for(i in runFolders){ dir.create(path = i, recursive = TRUE) }

    ####################
    # Build cube
    ####################
    # equal parameters between males and females
    #  all deposition is turned off
    #  no genotype-specific fitness costs
    cube <- MGDrivE2::cubeSEMcis(pF = paramCombo[[nP,'pF']], qF = paramCombo[[nP,'qF']], rF = paramCombo[[nP,'rF']],
                                 aF = paramCombo[[nP,'aF']], bF = paramCombo[[nP,'bF']], cF = paramCombo[[nP,'cF']],
                                 mmrF = paramCombo[[nP,'mmrF']])

    ####################
    # Setup Petri Net
    ####################
    # Places and transitions
    if(numPatch==1){
      # one node
      SPN_P <- MGDrivE2::spn_P_lifecycle_node(params = bioPars, cube = cube)
      SPN_T <- MGDrivE2::spn_T_lifecycle_node(spn_P = SPN_P, params = bioPars,
                                              cube = cube, feqTol = 1e-10)
    } else {
      # network
      SPN_P <- MGDrivE2::spn_P_lifecycle_network(num_nodes = numPatch, params = bioPars,
                                                 cube = cube)
      SPN_T <- MGDrivE2::spn_T_lifecycle_network(spn_P = SPN_P, params = bioPars,
                                                 cube = cube, m_move = adjMat,
                                                 feqTol = 1e-10)
    }

    # Stoichiometry matrix
    S <- MGDrivE2::spn_S(spn_P = SPN_P, spn_T = SPN_T)

    ####################
    # Equilibrium and Hazards
    ####################
    # calculate equilibrium and setup initial conditions
    initialCons <- MGDrivE2::equilibrium_lifeycle(params = bioPars, NF = numFem,
                                                  phi = bioPars$phi, log_dd = TRUE,
                                                  spn_P = SPN_P, cube = cube)

    # add movement stuff
    initialCons$params$mosquito_move_rates <- probRateList$gamma
    initialCons$params$mosquito_move_probs <- probRateList$mat

    # exact hazards for integer-valued state space
    exact_hazards <- MGDrivE2::spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube,
                                           params = initialCons$params, type = 'life',
                                           log = TRUE, exact = TRUE, tol = NaN,
                                           verbose = FALSE)

    ####################
    # Releases
    ####################
    # gene drive releases
    if(paramCombo[[nP, 'numRel1']] > 0){
      # setup que to release in first node, if there are multiple nodes
      relNode <- if( numPatch>1 ) '_1' else ''

      # make sure we're performing releases, otherwise skip everything
      gdEvents <- data.frame('var' = paste0('F_', cube$releaseType, '_', cube$wildType, relNode),
                             'time' = seq(from = relStart1,
                                          length.out = paramCombo[[nP, 'numRel1']],
                                          by = relInt1),
                             'value' = paramCombo[[nP, 'sizeRel1']],
                             'method' = 'add',
                             stringsAsFactors = FALSE)


      if(paramCombo[[nP,'numRel2']] > 0){
        # make sure we're doing SEM spraying, otherwise skip
        # sem "releases"
        # setup times for all of them
        timeVar <- seq.int(from = relStart2, length.out = paramCombo[[nP,'numRel2']], by = relInt2)

        # this is how the spn builds the statespace
        #  we have to recreate this perfectly, or the swapping won't work
        # females_unmated <- file.path("U",g, fsep = "_")
        # females_unmated <- file.path("U_", g, "_", node_id, fsep = "")
        #
        # males <- file.path("M",g,fsep = "_")
        # males <- file.path("M",g,node_id,fsep = "_")
        #
        # females <- file.path("F_",rep(g, each = nG),"_",g, fsep = "")
        # females <- file.path("F_",rep(g, each = nG),"_",g,"_",node_id, fsep = "")
        uEvents <- data.frame('var' = rep(file.path("U_", swapMat[ ,"off"], relNode, fsep = ''),
                                          each = paramCombo[[nP,'numRel2']]),
                              'value' = rep(file.path("U_", swapMat[ ,"on"], relNode, fsep = ''),
                                          each = paramCombo[[nP,'numRel2']]),
                              'time' = timeVar,
                              'method' = "swap")

        mEvents <- data.frame('var' = rep(file.path("M_", swapMat[ ,"off"], relNode, fsep = ''),
                                          each = paramCombo[[nP,'numRel2']]),
                              'value' = rep(file.path("M_", swapMat[ ,"on"], relNode, fsep = ''),
                                          each = paramCombo[[nP,'numRel2']]),
                              'time' = timeVar,
                              'method' = "swap")

        # grab all female states from petri net
        femHold <- grep(pattern = "^F_", x = SPN_P$u, value = TRUE)
        # if more than 1 node, we assume only node 1 gets sprayed
        #  this will fail if the nodeID string is duplicated in the genotypes somehow
        if(numPatch>1 ) femHold <- grep(pattern = relNode, x = femHold, fixed = TRUE, value = TRUE)
        # grap only the states with impacted genotypes
        femVar <- grep(pattern = paste0(swapMat[ ,"off"], collapse = "|"), x = femHold, value = TRUE)
        # swap appropriate genotypes
        femValue <- femVar
        for(i in 1:NROW(swapMat)){
          femValue <- gsub(pattern = swapMat[i,"off"], replacement = swapMat[i,"on"], x = femValue, fixed = TRUE)
        }
        # female events
        fEvents <- data.frame('var' = rep(femVar, each = paramCombo[[nP,'numRel2']]),
                              'value' = rep(femValue, each = paramCombo[[nP,'numRel2']]),
                              'time' = timeVar,
                              'method' = "swap")


        # combine all events for simulation
        events <- rbind(gdEvents, uEvents, mEvents, fEvents)

      } else {
        # no SEM releases, basically rename GD releases
        events <- gdEvents
      } # end SEM check

    } else {
      # we aren't doing GD releases, so no point in SEM spraying
      events <- NULL
    } # end GD release check


    ########################################
    ### Simulations
    ########################################
    ####################
    # Sims
    ####################
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
      parallel::clusterExport(cl = cl,
                              varlist = c('initialCons', 'simTime', 'sampTime',
                                          'dtStoch', 'runFolders', 'S',
                                          'exact_hazards', 'events')
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


    # run reps
    retList <- parallel::clusterApplyLB(cl = cl, x = 1:numRep, fun = function(x){
      MGDrivE2::sim_trajectory_CSV(
        x0 = initialCons$M0, tmax = simTime, dt = sampTime,
        dt_stoch = dtStoch, folders = runFolders[x], sampler = "tau",
        stage = c("M", "F"), S = S, Sout = NULL,
        hazards = exact_hazards, events = events, verbose = FALSE, maxhaz = 1e12
      )

      gc()
    })

    # stop cluster
    parallel::stopCluster(cl)

    ####################
    # Aggregate
    ####################
    # this is relatively fast, but if we're already running in parallel, lets make
    #  it faster
    if(Sys.info()[['sysname']] == "Windows"){
      # windows can't run fork clusters
      #  it has to use sockets
      cl <- parallel::makePSOCKcluster(names = ifelse(test = numCores>=2, yes = 2, no = numCores))

      # export required objects to each socket
      parallel::clusterExport(cl = cl,
                              varlist = c('runDir', 'SPN_P', 'simTime', 'sampTime'))

    } else {
      # *nix systems
      cl <- parallel::makeForkCluster(nnodes = ifelse(test = numCores>=2, yes = 2, no = numCores))
    }

    # split/aggregate
    # splits output by patch
    # aggregates females by mate genotype
    retList <- parallel::clusterApplyLB(cl = cl, x = c("M","FS"), fun = function(x){
      # this is very short
      MGDrivE2::split_aggregate_CSV(read_dir = runDir, spn_P = SPN_P,
                                    tmax = simTime, dt = sampTime, stage = x,
                                    erlang = FALSE, sum_fem = FALSE,
                                    verbose = TRUE, rem_file = TRUE)
    })

    # stop cluster
    parallel::stopCluster(cl)

    # print time it took
    print(paste0('  Param: ', capture.output(difftime(time1 = Sys.time(), time2 = startTimeParams))))

  } # end loop over parameter space

  # print time it took
  print(paste0(' Migration: ', capture.output(difftime(time1 = Sys.time(), time2 = startTimeLand))))

} # end loop over migration rates

# print time it took
print(paste0('All Sims: ', capture.output(difftime(time1 = Sys.time(), time2 = startTime))))


###############################################################################
### Analysis and Plotting
###############################################################################
########################################
### Analysis
########################################
# grab the parameter directories
migDirs <- list.dirs(path = baseOut, full.names = TRUE, recursive = FALSE)
migParamDirs <- lapply(X = migDirs, FUN = list.dirs, full.names = TRUE, recursive = FALSE)

# taken from the cube, so I didn't have to rerun it
genos <- c("WW","GW","UW","RW","VW","HW","SW","GG","GU","GR","GV","GH","GS","UU",
           "RU","UV","HU","SU","RR","RV","HR","RS","VV","HV","SV","HH","HS","SS")

# build the genotypes of interest
#  have to do this because there's no total by default, and I want to see that
goiList <- c(setNames(object = as.list(genos), nm = genos),
             list("Total"=genos))

# # example of counting alleles instead of genotypes
# goiList <- list("W"=c("WW","WW","GW","UW","RW","VW","HW","SW"),
#                 "G"=c("GW","GG","GG","GU","GR","GV","GH","GS"),
#                 "H"=c("HW","GH","HU","HR","HV","HH","HH","HS"),
#                 "U"=c("UW","GU","UU","UU","RU","UV","HU","SU"),
#                 "R"=c("RW","GR","RU","RR","RR","RV","HR","RS"),
#                 "V"=c("VW","GV","UV","RV","VV","VV","HV","SV"),
#                 "S"=c("SW","GS","SU","RS","SV","HS","SS","SS"),
#                 "Total"=c(genos, genos))

# loop over parameter directories, do analysis!
for(wDir in migParamDirs[[1]]){
  MGDrivE2::analyze_ggplot_CSV(read_dir = wDir, sex = "both", patch_agg = FALSE,
                               goi = goiList, drop_zero_goi = TRUE)
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
gPDat <- read.csv(file = file.path(migAnalysisDirs[[1]][6], "plot.csv"),
                  header = TRUE, sep = ",")
# columns: "Sex" "Patch" "Time" "GOI" "Min" "Mean" "Max" "2.5%" "50%" "97.5%"
gMDat <- read.csv(file = file.path(migAnalysisDirs[[1]][6], "metrics.csv"),
                  header = TRUE, sep = ",", check.names = FALSE)

####################
# Traces
####################
# data
ggplot(data = gPDat ) +
# color scheme and organization
geom_path(aes(x = Time, y = Count, group = interaction(Repetitions, GOI),
              color = GOI), alpha = 0.15) +
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

