################################################################################
#
#   MGDrivE2: Auxiliary data wrangling functions
#   Ndeffo Mbah Lab
#   Jared Bennett, jared_bennett@berkeley.edu
#   May 2022
#
################################################################################
# This has functions for analyzing adult data from V1 and V2, stored as .csv files
#


male/female = "m","f","both","agg"
read_dir <- "~/Desktop/OUTPUT/MGDrivE/V1/100/0900_0900_0050_0900_0900_0900_0900_0900_0050_00_1000"
read_dir <- "/home/gilchrist/Desktop/OUTPUT/MGDrivE/V1_3/100/0900_0900_0050_0900_0900_0900_0900_0900_0050_00_1000"


goi <- list("A" = c("WWWW","GWWW"), "B" = "SSEE")

jbFunc <- function(read_dir, write_dir = read_dir, sex = "both", patch_agg = FALSE,
                   goi = NULL, drop_zero_goi = TRUE){

  ##########
  # Checks
  ##########
  # required parameters
  if(missing(read_dir)){
    stop("Please provide: 'read_dir'")
  }

  # check read_dir
  if(!dir.exists(paths = read_dir)){
    stop("The 'read_dir' does not exist")
  }

  # check write_dir
  if(!dir.exists(paths = write_dir)){
    stop("The 'write_dir' does not exist, please create it")
  }

  # check sex
  mVal <- c('m','male')
  fVal <- c('f','female')
  if(!(sex %in% c(mVal,fVal,'both','agg'))){
    stop("'sex' must be one of: 'm', 'f', 'both', 'agg'")
  }
  # expand to more informative names
  #  let the default be "as-is"
  sexBool <- (sex == 'agg')
  sex <- switch(EXPR = sex, 'm'='male', 'f'='female',
                'both'=c('female','male'), 'agg'=c('female','male'))

  # check patch_agg
  if(!is.logical(patch_agg)){
    stop("'patch_agg' must be 'TRUE' or 'FALSE'")
  }

  # check drop_zer_gen
  if(!is.logical(drop_zero_gen)){
    stop("'drop_zero_gen' must be 'TRUE' or 'FALSE'")
  }

  ##########
  # Setup
  ##########
  # the search pattern subsets for male or female only input
  #  do this first in case someone has malformed data - then the appropriate sex
  #  only is checked.
  searchPattern <- c('^F_','^M_')
  if(length(sex) == 1){
    if(sex %in% mVal ){
      searchPattern <- '^M_'
    } else if(sex %in% fVal){
      searchPattern <- '^F_'
    }
  }

  # male/female file list
  # organization: repetitions -> male/female -> file name vec
  fileList <- lapply(X = list.dirs(path = read_dir, full.names = TRUE, recursive = FALSE),
                     FUN = function(x){lapply(X = searchPattern,
                                              FUN = list.files,
                                              path = x, full.names = TRUE)})

  # check input data
  if(length(fileList) == 0){
    stop("There are no directories in the 'read_dir', please check simulation output")
  }
  if(length(unique(unlist(lapply(X = fileList, FUN = lengths)))) != 1){
    stop("Not all of the repetitions have the same number of files, please check simulation output")
  }

  # read in example file
  testFile <- read.csv(file = fileList[[1]][[1]][1], header = TRUE)

  # get constants for later
  numReps <- length(fileList)
  numPatch <- length(fileList[[1]][[1]])
  genotypes <- colnames(testFile)[-1] # first column is "Time"
  numGeno <- length(genotypes)
  simTimes <- testFile[ ,"Time"]
  simTime <- length(simTimes)
  numSex <- length(sex)
  numCol <- numGeno+1L # because it reads in "Time", and then we drop it
  numRead <- simTime*numCol # simTime is fine because we skip 1 row, the header row

  # check goi and setup objects for later
  if(is.null(goi)) {
    # named list where each element is a genotype
    goi <- as.list(x = setNames(object = genotypes, nm = genotypes))
  } else if(!all(unlist(x = goi, use.names = FALSE) %in% genotypes)) {
    stop("elements of 'goi' must be valid genotypes in the files, please check goi list")
  }
  numGOI <- length(goi)
  goiNumeric <- lapply(X = goi, FUN = match, table = genotypes)
  goiLengths <- lengths(x = goi)


  ##########
  # data
  ##########
  # setup generic data array
  #  this may get reduced later by sex and/or patch
  dataArray <- array(data = 0, dim = c(numReps, numSex, numPatch, simTime, numGOI),
                     dimnames = list("Repetitions"=seq_len(numReps),
                                     "Sex"=sex,
                                     "Patch"=seq_len(numPatch),
                                     "Time"=simTimes,
                                     "GOI"=names(goi))
  )

  # wR = which rep
  # wS = which sex
  # wP = which patch
  # wGOIN = which goiNumeric
  for(wR in seq_len(numReps)){
    for(wS in seq_len(numSex)){
      for(wP in seq_len(numPatch)){
        # read in file, drop header
        holdMat <- matrix(data = scan(file = fileList[[wR]][[wS]][[wP]],
                                      what = numeric(), n = numRead,
                                      sep = ",", skip = 1, quiet = TRUE),
                          nrow = simTime, ncol = numCol, byrow = TRUE)[ ,-1]

        # collapse genotypes into genotypes of interest
        for(wGOIN in seq_len(numGOI)){
          dataArray[wR,wS,wP, ,wGOIN] <- .rowSums(x = holdMat[ ,goiNumeric[[wGOIN]],drop=FALSE],
                                                  m = simTime, n = goiLengths[wGOIN])
        }

      } # end patch loop
    } # end sex loop
  } # end rep loop

  # handle patch reduction if applicable
  if(patch_agg){
    # is there a better way to do this?
    #  it seems really bad
    # turn array into a list of arrays, slicing along the patch dimension
    # aggregate the list back into a single array
    dataArray <- Reduce(f = "+",
                        x = lapply(X = seq_len(dim(dataArray)[3]),
                                   FUN = function(x)dataArray[ , ,x, , ,drop = FALSE]))
    # update dimension name
    dimnames(dataArray)$Patch <- "all"
  } # end patch reduction

  # handle sex aggregation if applicable
  if(sexBool){
    # combine over sex
    dataArray <- dataArray[ ,1, , , ,drop = FALSE] + dataArray[ ,2, , , ,drop = FALSE]
    # update dimension name
    dimnames(dataArray)$Sex <- "agg"
  } # end sex reduction

  # handle GOI reduction if applicable
  if(drop_zero_goi){
    # get indices of non-zero GOI
    keepIdx <- which(vapply(X = seq_len(numGOI),
                            FUN = function(x){sum(dataArray[ , , , ,x])>=1},
                            FUN.VALUE = logical(length = 1)))
    # safety check
    if(length(keepIdx)==0){
      stop("There are no genotypes of interest in the analysis, please check simulation output")
    }
    # subset, handles names internally
    dataArray <- dataArray[ , , , ,keepIdx,drop = FALSE]
  } # end drop check

















}

































################################################################################
# Main Split-Aggregate Function
################################################################################

#' Split CSV output by Patch and Aggregate by Mate or Dwell-Stage
#'
#' This function reads in the output files from \code{\link{sim_trajectory_CSV}}
#' and splits them into smaller files. The files are output by patch, with the
#' appropriate patch numbers for mosquitoes or humans, and specific stages are
#' aggregated by a given metric. \cr
#'
#' Given the \code{read_dir}, this function assumes the follow file structure: \cr
#'  * read_dir
#'    * repetition 1
#'      * M.csv
#'      * FS.csv
#'      * ... \cr
#'    * repetition 2
#'      * M.csv
#'      * FS.csv
#'      * ... \cr
#'    * repetition 3
#'    * ... \cr
#'
#' This function expects the \code{write_dir} to be empty, and it sets up the
#' same file structure as the \code{read_dir}. For a 2-node simulation, the output
#' will be organized similar to: \cr
#'  * write_dir
#'    * repetition 1
#'      * M_0001.csv
#'      * M_0002.csv
#'      * FS_0001.csv
#'      * FS_0001.csv
#'      * ... \cr
#'    * repetition 2
#'      * M_0001.csv
#'      * M_0002.csv
#'      * FS_0001.csv
#'      * FS_0001.csv
#'      * ... \cr
#'    * repetition 3
#'    * ... \cr
#'
#' \code{stage} defines which life-stages the function will analyze. These stages
#' must be any combination of: "E", "L", "P", "M", "U", "FS", "FE", "FI", "H".
#' These must come from the set of stages provided to \code{\link{sim_trajectory_CSV}}
#' via the \code{stage} argument. It can be less than what was printed by the simulation,
#' but any extra stages provided, but not printed, will throw a warning and then
#' be ignored.
#'
#' \code{erlang} defines how aquatic (eggs, larvae, and pupae) stages and adult females
#' (only mated females) are aggregated. By default, \code{erlang} is FALSE, and
#' all of these stages are summarized by genotype only, combining any Erlang-distributed
#' dwell stages (for eggs, larvae, and pupae) or latent infection (for adult females)
#' stages. If \code{erlang} is TRUE, summaries are returned by dwell stage or infection
#' status, combining any genotype information. \cr
#' Female summaries always combine over mate-genotype, so only female genotypes
#' are returned.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' \code{tmax}, \code{dt} define the last sampling
#' time, and each sampling time in-between.
#'
#' For more details about using this function to process CSV output see:
#' \code{vignette("data-analysis", package = "MGDrivE2")}
#'
#' @param read_dir Directory where output was written to
#' @param write_dir Directory to write output to. Default is read_dir
#' @param stage Life stage to print, see details
#' @param spn_P Places object, see details
#' @param tmax The final time to end simulation
#' @param dt The time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param erlang Boolean, default is FALSE, to return summaries by genotype
#' @param sum_fem if \code{TRUE}, in addition to FS, FE, FI output by node and repetition, output an
#'                additional file F which sums over infection states (S,E,I). Does nothing if the
#'                simulation did not include epi dynamics.
#' @param rem_file Remove original output? Default is FALSE
#' @param verbose Chatty? Default is TRUE
#'
#' @return Writes output to files in write_dir
#'
#' @importFrom utils write.table
#'
#' @export
split_aggregate_CSV <- function(
  read_dir,
  write_dir = read_dir,
  stage = c("E","L","P","M","U","FS","FE","FI","H"),
  spn_P, tmax, dt,
  erlang = FALSE,
  sum_fem = FALSE,
  rem_file=FALSE,
  verbose=TRUE
){

  ##########
  # Checks
  ##########
  # required parameters
  if(any(c(missing(read_dir),missing(spn_P),missing(tmax),missing(dt)))){
    stop("Please provide 'read_dir', 'spn_P', 'tmax', and 'dt'.")
  }
  t0 <- 0
  tt <- tmax

  # check read_dir
  if(!dir.exists(read_dir)){
    stop("The 'read_dir' does not exist.")
  }

  # check write_dir
  if(!dir.exists(paths = write_dir)){
    stop("The 'write_dir' does not exist, please create it.")
  }

  # check stage
  if(!all(stage %in% c("E","L","P","M","U","FS","FI","FE","H"))){
    stop("Analysis stages must match files output. ",
         "Please provide a combination of: 'E','L','P','M','U','FS','FI','FE','H' ")
  }
  stage <- unique(x = stage)


  ##########
  # Input and Output Setup
  ##########
  # setup output directories
  repDirs <- file.path(write_dir, list.dirs(path = read_dir, recursive = FALSE, full.names = FALSE))
  for(wDir in repDirs){dir.create(path = wDir, showWarnings = FALSE, recursive = FALSE)}

  # get all files to work on, sorted by type
  fileList <- lapply(X = stage, FUN = function(x){
                      list.files(path = read_dir, pattern = paste0("^",x),
                                 full.names = TRUE, recursive = TRUE)
                     })

  # check for things that aren't there, remove them
  zeroIdx <- which(!lengths(x = fileList))
  if(length(zeroIdx)>0){
    # warning
    warning(stage[zeroIdx], " was/were not output by the simulation, ",
            "and will be removed from the analysis.\n")
    # remove things
    fileList[zeroIdx] <- NULL
    stage <- stage[-zeroIdx]
  }


  ##########
  # Derived Parameters
  ##########
  # get nodes
  #  Since some simulations have humans, mosquitoes, and mixed nodes, use the
  #  places object to determine those numbers.
  #  use of "egg" was a choice, any of the mosquito life stages would be fine
  mosyNodes <- which(x = vapply(X = spn_P$ix,
                                FUN = function(x){!is.null(x$egg)},
                                FUN.VALUE = logical(length = 1)) )
  numMosyNodes <- length(mosyNodes)

  hNodes <- which(x = vapply(X = spn_P$ix,
                                FUN = function(x){!is.null(x$humans)},
                                FUN.VALUE = logical(length = 1)) )
  numHNodes <- length(hNodes)

  # get genotype info
  mosyGenos <- colnames(spn_P$ix[[mosyNodes[1]]]$egg)
  numMosyGenos <- length(mosyGenos)
  hGenos <- switch((length(spn_P$ix[[hNodes[1]]]$humans) == 2)+1, c("S","E","I","R"), c("S","I"))
  numHGenos <- length(hGenos)

  # get length of Erlang stages
  nELP <- c("E" = NROW(spn_P$ix[[mosyNodes[1]]]$egg),
            "L" = NROW(spn_P$ix[[mosyNodes[1]]]$larvae),
            "P" = NROW(spn_P$ix[[mosyNodes[1]]]$pupae) )
  # numEIP is female latent stages, if there isn't infection, there is only 1 "stage"
  #  if there is infection, there is >=3 states, S, E...., I, so numEIP subtract 2 for S and  I
  numEIP <- ifelse(dim(spn_P$ix[[mosyNodes[1]]]$females)[3] == 1, 1, dim(spn_P$ix[[mosyNodes[1]]]$females)[3]-2)

  # if there's no epi dynamics, rename FS -> F
  if(numEIP == 1){
    stage[which(stage == "FS")] <- "F"
  }

  # sampling times, and number of samples in files (ie, number of rows)
  times <- seq(from=t0,to=tt,by=dt)
  nTimes <- length(times)

  numStage <- length(stage)

  # node names for output
  mosyNodeNames <- formatC(x = mosyNodes, width = 4, format = "d", flag = "0")
  hNodeNames <- formatC(x = hNodes, width = 4, format = "d", flag = "0")


  ##########
  # Initialize Text
  ##########
  if(verbose){
    # what we're doing
    cat("Analyzing", numStage, "directories.\n")

    # whether or not we remove files
    if(rem_file){
      cat("\tRemoving original files.\n\n")
    } else {
      cat("\tNot removing original files.\n\n")
    }
  }


  ##########
  # Dispatch Loop
  ##########
  for(analysis in 1:numStage){

    # generic file names that work for most things here
    #  need to redefine for humans and erlang female summaries
    outputNames <- lapply(X = repDirs, FUN = function(x){
                          file.path(x,file.path(stage[analysis],"_",mosyNodeNames,".csv",fsep = ""))})

    # indication that something is working
    if(verbose){
      cat("Starting analysis", analysis, "of", numStage, " ...  ")
    }


    # begin logic tree
    if(stage[analysis] %in% c("E","L","P")){
      # check for erlang or geno summary
      if(erlang){
        base_erlang(fileVec = fileList[[analysis]],
                    outList = lapply(X = repDirs, FUN = function(x){file.path(x,file.path(stage[analysis],"E_",mosyNodeNames,".csv",fsep = ""))}),
                    genos = mosyGenos, nGenos = numMosyGenos, nErlang = nELP[stage[analysis]],
                    times = times, nTimes = nTimes, nNodes = numMosyNodes)
      } else {
        base_gen(fileVec = fileList[[analysis]], outList = outputNames,
                 genos = mosyGenos,nGenos = numMosyGenos,nIDX1 = nELP[stage[analysis]],
                 times = times,nTimes = nTimes,nNodes = numMosyNodes)
      }

    } else if(stage[analysis] %in% c("M","U")){
      base_MUH(fileVec = fileList[[analysis]], outList = outputNames, genos = mosyGenos,
               nGenos = numMosyGenos, nTimes = nTimes, nNodes = numMosyNodes)

    } else if(stage[analysis] %in% c("FS","FI","F")){
      # check for erlang or geno summary
      if(erlang){
        if(verbose){
          cat("Skip\n")
        }
        next
      } else {
        base_gen(fileVec = fileList[[analysis]], outList = outputNames,
                 genos = mosyGenos, nGenos = numMosyGenos, nIDX1 = numMosyGenos,
                 times = times, nTimes = nTimes, nNodes = numMosyNodes)
      }

    } else if(stage[analysis] == "FE"){
      if(erlang){
        # need all 3 female files, so get index for them
        fIDX <- match(x = c("FS","FE","FI"), table = stage)
        if(any(is.na(fIDX))){
          warning("\nErlang summary of adult females requires 'FS', 'FE', and 'FI' files.\n",
                  c("FS","FE","FI")[is.na(fIDX)], " files were not found.\n",
                  "\t Skipping Erlang summary of adult females.")
          next
        }

        base_erlang_F(fileList = fileList[fIDX],
                      outList = lapply(X = repDirs, FUN = function(x){file.path(x,file.path("FSEI_",mosyNodeNames,".csv",fsep = ""))}),
                      nGenos = numMosyGenos, nErlang = numEIP, times = times,
                      nTimes = nTimes, nNodes = numMosyNodes)

      } else {
        base_gen_FE(fileVec = fileList[[analysis]], outList = outputNames,
                    genos = mosyGenos, nGenos = numMosyGenos, nIDX1 = numEIP,
                    times = times, nTimes = nTimes, nNodes = numMosyNodes)
      }

    } else if(stage[analysis] == "H"){
      # have to redefine filenames for humans
      base_MUH(fileVec = fileList[[analysis]],
               outList = lapply(X = repDirs, FUN = function(x){file.path(x,file.path("H_",hNodeNames,".csv",fsep = ""))}),
               genos = hGenos, nGenos = numHGenos, nTimes = nTimes, nNodes = numHNodes)
    } # end if tree
  } # end analysis loop

  # check if we need to sum all females (SEI)
  if(sum_fem & (numEIP != 1)){
    fIDX <- match(x = c("FS","FE","FI"), table = stage)
    if(any(is.na(fIDX))){
      warning("\nSumming of adult females requires 'FS', 'FE', and 'FI' files.\n",
              c("FS","FE","FI")[is.na(fIDX)], " files were not found.\n",
              "\t Skipping summary of adult females.")
    } else {
      base_sum_F(fileList = fileList[fIDX],
                    outList = lapply(X = repDirs, FUN = function(x){file.path(x,file.path("F_",mosyNodeNames,".csv",fsep = ""))}),
                    genos = mosyGenos,nGenos = numMosyGenos, nErlang = numEIP, times = times,
                    nTimes = nTimes, nNodes = numMosyNodes)
    }
  }

  # done with analysis
  if(verbose){
    cat("Done\n")
  }

  # check if removing original output
  if(rem_file){file.remove(unlist(fileList))}

} # end function


################################################################################
# Base Split-Aggregate Functions
################################################################################
#######################################
# Genotype Only Stages
#######################################

#' Base Summary for Males, Unmated Females, and Humans
#'
#' This function takes a given stage (males, unmated females, or humans) and
#' summarizes them by genotype (infection status for humans), writing output to
#' provided folders.
#'
#' This function is a base function used in \code{\link{split_aggregate_CSV}}.
#'
#' @param fileVec Vector of files for analysis
#' @param outList List of files, organized by repetition, to write output
#' @param genos Genotypes to summarize by
#' @param nGenos Number of genotypes
#' @param nTimes Number of sampled times
#' @param nNodes Number of nodes in the network
#'
#' @return None
#'
base_MUH <- function(fileVec,outList,genos,nGenos,nTimes,nNodes){

  # indexing
  #  each file gets nGenos of stuff, plus the "Time" column, so 1:nGenos.
  #    since "Time" is always the first column, shift the rest by 1
  #  Each further node must step over nGeno number of previous columns,
  #  so steps by nGenos, but the first set doesn't step, shift 0:(nNodes-1)
  gen_node_idx <- outer(X = 1:nGenos + 1, Y = 0:(nNodes-1)*nGenos, FUN = "+")

  # column names for output
  columnNames <- c("Time",genos)

  # loop variables
  numCol <- nGenos*nNodes + 1
  numTot <- numCol*nTimes

  # loop over repetitions
  for(r in 1:length(fileVec)){

    # read in files
    fileIn = matrix(data = scan(file = fileVec[r], what = numeric(), sep = ",",
                              n = numTot, skip = 1, quiet = TRUE),
                    ncol = numCol, nrow = nTimes, byrow = TRUE)

    # loop over nodes and write output
    for(node in 1:nNodes){

      # subset and write output
      write.table(x = fileIn[ ,c(1,gen_node_idx[ ,node])],
                  file = outList[[r]][node], sep = ",", row.names = FALSE,
                  col.names = columnNames, quote = FALSE)

    } # end loop over nodes
  } # end loop over reps

} # end function


#######################################
# Genotype and 1 Other Index
#######################################

#' Base Summary for Eggs, Larvae, Pupae, Susceptible Females, and Infectious Females
#'
#' This function takes a given stage and summarizes them by genotype, writing
#' output to provided folders.
#'
#' This function is a base function used in \code{\link{split_aggregate_CSV}}.
#'
#' @param fileVec Vector of files for analysis
#' @param outList List of files, organized by repetition, to write output
#' @param genos Genotypes to summarize by
#' @param nGenos Number of genotypes
#' @param nIDX1 First index to expand over, nE/nL/nP for aquatic stages, 1 for the rest
#' @param times Vector of sampling times
#' @param nTimes Number of sampled times
#' @param nNodes Number of nodes in the network
#'
#' @return None
#'
base_gen <- function(fileVec,outList,genos,nGenos,nIDX1,times,nTimes,nNodes){

  # indexing
  #  array, indexes by first index, then genotype, then node
  #  add 1 to skip "Time" index
  aggIDX <- outer(X = outer(X = 1:nIDX1+1, Y = 0:(nGenos-1)*nIDX1, FUN = "+"),
                  Y = 0:(nNodes-1)*nGenos*nIDX1, FUN = "+")

  # setup output
  retMat <- matrix(data = 0,nrow = nTimes, ncol = nGenos+1,
                   dimnames = list(NULL,c("Time",genos)) )
  retMat[ ,1] <- times

  # loop variables
  numCol <- nIDX1*nGenos*nNodes + 1
  numTot <- numCol*nTimes

  # loop over repetitions
  for(r in 1:length(fileVec)){

    # read in files
    fileIn = matrix(data = scan(file = fileVec[r], what = numeric(), sep = ",",
                              n = numTot, skip = 1, quiet = TRUE),
                    ncol = numCol, nrow = nTimes, byrow = TRUE)

    # loop over nodes and write output
    for(node in 1:nNodes){
      # loop over aggregation
      for(agg in 1:nGenos){
        retMat[ ,agg+1] <- rowSums(x = fileIn[ ,aggIDX[ ,agg,node], drop = FALSE])
      }

      # write output
      write.table(x = retMat, file = outList[[r]][node], sep = ",",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)

    } # end loop over nodes
  } # end loop over reps

} # end function


#######################################
# Genotype - Latent Females
#######################################

#' Base Summary for Latent Females
#'
#' This function takes 'E' stage females and summarizes them by genotype, writing
#' output to provided folders.
#'
#' This function is a base function used in \code{\link{split_aggregate_CSV}}.
#'
#' @param fileVec Vector of files for analysis
#' @param outList List of files, organized by repetition, to write output
#' @param genos Genotypes to summarize by
#' @param nGenos Number of genotypes
#' @param nIDX1 First index to expand over, nE/nL/nP for aquatic stages, 1 for the rest
#' @param times Vector of sampling times
#' @param nTimes Number of sampled times
#' @param nNodes Number of nodes in the network
#'
#' @return None
#'
base_gen_FE <- function(fileVec,outList,genos,nGenos,nIDX1,times,nTimes,nNodes){

  # indexing
  #  array, indexes by first genotype, then genotype, then numEIP, then node
  #  add 1 to skip "Time" index
  aggIDX <- outer(X = outer(X = outer(X = 1:nGenos+1, Y = 0:(nGenos-1)*nGenos, FUN = "+"),
                            Y = 0:(nIDX1-1)*nGenos*nGenos, FUN = "+"),
                  Y = 0:(nNodes-1)*nGenos*nGenos*nIDX1, FUN = "+")

  # setup output
  retMat <- matrix(data = 0,nrow = nTimes, ncol = nGenos+1,
                   dimnames = list(NULL,c("Time",genos)) )
  retMat[ ,1] <- times

  # loop variables
  numCol <- nGenos*nGenos*nIDX1*nNodes + 1
  numTot <- numCol*nTimes

  # loop over repetitions
  for(r in 1:length(fileVec)){

    # read in files
    fileIn <- matrix(data = scan(file = fileVec[r], what = numeric(), sep = ",",
                              n = numTot, skip = 1, quiet = TRUE),
                    ncol = numCol, nrow = nTimes, byrow = TRUE)

    # loop over nodes and write output
    for(node in 1:nNodes){
      # loop over aggregation
      for(agg in 1:nGenos){
        retMat[ ,agg+1] <- rowSums(x = fileIn[ ,aggIDX[ ,agg, ,node], drop = FALSE])
      }

      # write output
      write.table(x = retMat, file = outList[[r]][node], sep = ",",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)

    } # end loop over nodes
  } # end loop over reps

} # end function


#######################################
# Erlang - Aquatic
#######################################

#' Base Summary of Erlang Stages for Aquatic Life Stages
#'
#' This function takes the given aquatic stage and summarizes them by Erlang-distributed
#' dwell times, writing output to provided folders.
#'
#' This function is a base function used in \code{\link{split_aggregate_CSV}}.
#'
#' @param fileVec Vector of files for analysis
#' @param outList List of files, organized by repetition, to write output
#' @param genos Genotypes to summarize by
#' @param nGenos Number of genotypes
#' @param nErlang Number of Erlang stages
#' @param times Vector of sampling times
#' @param nTimes Number of sampled times
#' @param nNodes Number of nodes in the network
#'
#' @return None
#'
base_erlang <- function(fileVec,outList,genos,nGenos,nErlang,times,nTimes,nNodes){

  # indexing
  #  array, indexes by first index, then genotype, then node
  #  add 1 to skip "Time" index
  aggIDX <- outer(X = outer(X = 1:nErlang+1, Y = 0:(nGenos-1)*nErlang, FUN = "+"),
                  Y = 0:(nNodes-1)*nGenos*nErlang, FUN = "+")

  # setup output
  retMat <- matrix(data = 0,nrow = nTimes, ncol = nErlang+1,
                   dimnames = list(NULL,c("Time",1:nErlang)) )
  retMat[ ,1] <- times

  # loop variables
  numCol <- nErlang*nGenos*nNodes + 1
  numTot <- numCol*nTimes

  # loop over repetitions
  for(r in 1:length(fileVec)){

    # read in files
    fileIn <- matrix(data = scan(file = fileVec[r], what = numeric(), sep = ",",
                              n = numTot, skip = 1, quiet = TRUE),
                    ncol = numCol, nrow = nTimes, byrow = TRUE)

    # loop over nodes and write output
    for(node in 1:nNodes){
      # loop over aggregation
      for(agg in 1:nErlang){
        retMat[ ,agg+1] <- rowSums(x = fileIn[ ,aggIDX[agg, ,node], drop = FALSE])
      }

      # write output
      write.table(x = retMat, file = outList[[r]][node], sep = ",",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)

    } # end loop over nodes
  } # end loop over reps

} # end function


#######################################
# Erlang - Adult Females
#######################################

#' Base Summary of Erlang Stages for Adult Females
#'
#' This function takes ALL of the adult female stages and summarized them by
#' Erlang-distributed latent infection, writing output to provided folders.
#'
#' This function is a base function used in \code{\link{split_aggregate_CSV}}.
#'
#' @param fileList Length 3 list holding 'FS','FE', and 'FI' files for analysis
#' @param outList List of files, organized by repetition, to write output
#' @param nGenos Number of genotypes
#' @param nErlang Number of Erlang stages
#' @param times Vector of sampling times
#' @param nTimes Number of sampled times
#' @param nNodes Number of nodes in the network
#'
#' @return None
#'
base_erlang_F <- function(fileList,outList,nGenos,nErlang,times,nTimes,nNodes){

  # indexing
  #  array, indexes by first genotype, then genotype, then numEIP, then node
  #  add 1 to skip "Time" index
  aggE <- outer(X = outer(X = outer(X = 1:nGenos+1, Y = 0:(nGenos-1)*nGenos, FUN = "+"),
                          Y = 0:(nErlang-1)*nGenos*nGenos, FUN = "+"),
                Y = 0:(nNodes-1)*nGenos*nGenos*nErlang, FUN = "+")

  aggSI <- outer(X = outer(X = 1:nGenos+1, Y = 0:(nGenos-1)*nGenos, FUN = "+"),
                 Y = 0:(nNodes-1)*nGenos*nGenos, FUN = "+")


  # setup output
  #  ncol = Time, S, E(1:numEIP), I
  retMat <- matrix(data = 0,nrow = nTimes, ncol = nErlang+3,
                   dimnames = list(NULL,c("Time","S",paste0("E",1:nErlang),"I")) )
  retMat[ ,1] <- times

  # loop variables
  numColSI <- nGenos*nGenos*nNodes + 1
  numTotSI <- numColSI*nTimes
  numColE <- nGenos*nGenos*nErlang*nNodes + 1
  numTotE <- numColE*nTimes
  iIndex <- nErlang + 3

  # loop over repetitions
  for(r in 1:length(fileList[[1]])){

    # read in files
    sFile <- matrix(data = scan(file = fileList[[1]][r], what = numeric(), sep = ",",
                                n = numTotSI, skip = 1, quiet = TRUE),
                    ncol = numColSI, nrow = nTimes, byrow = TRUE)

    eFile <- matrix(data = scan(file = fileList[[2]][r], what = numeric(), sep = ",",
                                n = numTotE, skip = 1, quiet = TRUE),
                    ncol = numColE, nrow = nTimes, byrow = TRUE)

    iFile <- matrix(data = scan(file = fileList[[3]][r], what = numeric(), sep = ",",
                                n = numTotSI, skip = 1, quiet = TRUE),
                    ncol = numColSI, nrow = nTimes, byrow = TRUE)


    # loop over nodes and write output
    for(node in 1:nNodes){
      # Aggregate by Erlang Stage
      # S
      retMat[ ,2] <- rowSums(x = sFile[ ,aggSI[ , ,node], drop = FALSE])

      # E
      for(agg in 1:nErlang){
        retMat[ ,agg+2] <- rowSums(x = eFile[ ,aggE[ , ,agg,node], drop = FALSE])
      }

      # I
      retMat[ ,iIndex] <- rowSums(x = iFile[ ,aggSI[ , ,node], drop = FALSE])


      # write output
      write.table(x = retMat, file = outList[[r]][node], sep = ",",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)

    } # end loop over nodes
  } # end loop over reps

} # end function


#######################################
# Sum all Adult Females (S,E,I) -> F
#######################################

#' Base Summary of Infection (SEI) Stages for Adult Females
#'
#' This function takes ALL of the adult female stages and summarized them by
#' Erlang-distributed latent infection, writing output to provided folders.
#'
#' This function is a base function used in \code{\link{split_aggregate_CSV}}.
#'
#' @param fileList Length 3 list holding 'FS','FE', and 'FI' files for analysis
#' @param outList List of files, organized by repetition, to write output
#' @param genos Genotypes to summarize by
#' @param nGenos Number of genotypes
#' @param nErlang Number of Erlang stages
#' @param times Vector of sampling times
#' @param nTimes Number of sampled times
#' @param nNodes Number of nodes in the network
#'
#' @return None
base_sum_F <- function(fileList,outList,genos,nGenos,nErlang,times,nTimes,nNodes){

  # indexing
  #  array, indexes by first genotype, then genotype, then numEIP, then node
  #  add 1 to skip "Time" index
  aggE <- outer(X = outer(X = outer(X = 1:nGenos+1, Y = 0:(nGenos-1)*nGenos, FUN = "+"),
                            Y = 0:(nErlang-1)*nGenos*nGenos, FUN = "+"),
                  Y = 0:(nNodes-1)*nGenos*nGenos*nErlang, FUN = "+")

  aggSI <- outer(X = outer(X = 1:nGenos+1, Y = 0:(nGenos-1)*nGenos, FUN = "+"),
                 Y = 0:(nNodes-1)*nGenos*nGenos, FUN = "+")


   # setup output
   retMat <- matrix(data = 0,nrow = nTimes, ncol = nGenos+1,
                    dimnames = list(NULL,c("Time",genos)) )
   retMat[ ,1] <- times


  # loop variables
  numColSI <- nGenos*nGenos*nNodes + 1
  numTotSI <- numColSI*nTimes
  numColE <- nGenos*nGenos*nErlang*nNodes + 1
  numTotE <- numColE*nTimes

  # loop over repetitions
  for(r in 1:length(fileList[[1]])){

    # read in files
    sFile <- matrix(
      data = scan(file = fileList[[1]][r], what = numeric(), sep = ",",
      n = numTotSI, skip = 1, quiet = TRUE),
      ncol = numColSI, nrow = nTimes, byrow = TRUE
    )

    eFile <- matrix(
      data = scan(file = fileList[[2]][r], what = numeric(), sep = ",",
      n = numTotE, skip = 1, quiet = TRUE),
      ncol = numColE, nrow = nTimes, byrow = TRUE
    )

    iFile <- matrix(
      data = scan(file = fileList[[3]][r], what = numeric(), sep = ",",
      n = numTotSI, skip = 1, quiet = TRUE),
      ncol = numColSI, nrow = nTimes, byrow = TRUE
    )

    # clear output
    retMat[ ,-1] <- retMat[ ,-1] * 0

    # loop over nodes and write output
    for(node in 1:nNodes){
      # loop over aggregation
      for(agg in 1:nGenos){
        retMat[ ,agg+1] <- retMat[ ,agg+1] + rowSums(x = sFile[ ,aggSI[ ,agg,node], drop = FALSE])
        retMat[ ,agg+1] <- retMat[ ,agg+1] + rowSums(x = eFile[ ,aggE[ ,agg, ,node], drop = FALSE])
        retMat[ ,agg+1] <- retMat[ ,agg+1] + rowSums(x = iFile[ ,aggSI[ ,agg,node], drop = FALSE])
      }

      # write output
      write.table(x = retMat, file = outList[[r]][node], sep = ",",
                  row.names = FALSE, col.names = TRUE, quote = FALSE)

    } # end loop over nodes


  } # end loop over reps

} # end function


################################################################################
# Main Statistical Summary Function
################################################################################

#' Summary Statistics for MGDrivE2
#'
#' This function reads in all repetitions for each patch and calculates either
#' the mean, quantiles, or both. User chooses the quantiles, up to 4 decimal places,
#' and enters them as a vector. Quantiles are calculated empirically. (order does not matter)  \cr
#'
#' Given the read_dir, this function assumes the follow file structure: \cr
#'  * read_dir
#'    * repetition 1
#'      * M_0001.csv
#'      * M_0002.csv
#'      * FS_0001.csv
#'      * FS_0001.csv
#'      * ... \cr
#'    * repetition 2
#'      * M_0001.csv
#'      * M_0002.csv
#'      * FS_0001.csv
#'      * FS_0001.csv
#'      * ... \cr
#'    * repetition 3
#'    * ... \cr
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' \code{t0}, \code{tt}, \code{dt} define the first sampling time, the last sampling
#' time, and each sampling time in-between.
#'
#' Output files are *.csv and contain the mean or quantile in the file name, e.g.
#' {stage}_Mean_(patchNum).csv and {stage}_Quantile_(quantNum)_(patchNum).csv.
#'
#' For more details about using this function to process CSV output see:
#' \code{vignette("data-analysis", package = "MGDrivE2")}
#'
#'
#' @param read_dir Directory to find repetition folders in
#' @param write_dir Directory to write output
#' @param mean Boolean, calculate mean or not. Default is TRUE
#' @param quantiles Vector of quantiles to calculate. Default is NULL
#' @param spn_P Places object, see details
#' @param tmax The final time to end simulation
#' @param dt The time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param rem_file Remove original output? Default is FALSE
#' @param verbose Chatty? Default is TRUE
#'
#' @return Writes output to files in write_dir
#'
#' @export
summarize_stats_CSV <- function(
  read_dir, write_dir=read_dir, mean=TRUE, quantiles=NULL,
  spn_P, tmax, dt, rem_file=FALSE, verbose=TRUE
){

  ##########
  # Checks
  ##########
  # required parameters
  if(any(c(missing(read_dir),missing(spn_P),missing(tmax),missing(dt)))){
    stop("Please provide 'read_dir', 'spn_P', 'tmax', and 'dt'.")
  }
  t0 <- 0
  tt <- tmax

  # check read_dir
  if(!dir.exists(read_dir)){
    stop("The 'read_dir' does not exist.")
  }

  # check write_dir
  if(!dir.exists(paths = write_dir)){
    stop("The 'write_dir' does not exist, please create it.")
  }

  # what analysis to do
  if(!mean && is.null(quantiles)){
    stop("User needs to specify the mean or which quantiles to calculate. ")
  }


  ##########
  # Input and Output Setup
  ##########
  # get all files to work on, sorted by type
  stage = c("E","L","P","M","U","FS","FE","FI","H","EE","LE","PE","FSEI","F")
  repDirs <- list.dirs(path = read_dir, recursive = FALSE, full.names = TRUE)
  fileList <- lapply(X = stage, FUN = function(x){
                                lapply(X = repDirs, FUN = list.files, pattern = paste0("^",x,"_"),
                                       full.names = TRUE, recursive = FALSE)
                     })

  # check for things that aren't there, remove them
  zeroIdx <- which(!vapply(X = lapply(X = fileList, FUN = lengths),
                           FUN = sum, FUN.VALUE = numeric(length = 1)))
  fileList[zeroIdx] <- NULL
  stage <- stage[-zeroIdx]


  ##########
  # Derived Parameters
  ##########
  # get nodes
  #  Since some simulations have humans, mosquitoes, and mixed nodes, use the
  #  places object to determine those numbers.
  #  use of "egg" was a choice, any of the mosquito life stages would be fine
  mosyNodes <- which(x = vapply(X = spn_P$ix,
                                FUN = function(x){!is.null(x$egg)},
                                FUN.VALUE = logical(length = 1)) )
  numMosyNodes <- length(mosyNodes)

  hNodes <- which(x = vapply(X = spn_P$ix,
                                FUN = function(x){!is.null(x$humans)},
                                FUN.VALUE = logical(length = 1)) )
  numHNodes <- length(hNodes)

  # get genotype info
  mosyGenos <- colnames(spn_P$ix[[mosyNodes[1]]]$egg)
  numMosyGenos <- length(mosyGenos)
  hGenos <- switch((length(spn_P$ix[[hNodes[1]]]$humans) == 2)+1, c("S","E","I","R"), c("S","I"))
  numHGenos <- length(hGenos)

  # get length of Erlang stages
  nELP <- c("E" = NROW(spn_P$ix[[mosyNodes[1]]]$egg),
            "L" = NROW(spn_P$ix[[mosyNodes[1]]]$larvae),
            "P" = NROW(spn_P$ix[[mosyNodes[1]]]$pupae) )
  # nEIP is female latent stages, if there isn't infection, there is only 1 "stage"
  #  if there is infection, there is >=3 states, S, E...., I, so nEIP subtract 2 for S and  I
  nEIP <- ifelse(dim(spn_P$ix[[mosyNodes[1]]]$females)[3] == 1, 1, dim(spn_P$ix[[mosyNodes[1]]]$females)[3]-2)

  # sampling times, and number of samples in files (ie, number of rows)
  times <- seq(from=t0,to=tt,by=dt)
  nTimes <- length(times)

  # node names for output
  mosyNodeNames <- formatC(x = mosyNodes, width = 4, format = "d", flag = "0")
  hNodeNames <- formatC(x = hNodes, width = 4, format = "d", flag = "0")


  # derived parameters specific to this function
  numStage <- length(stage)
  numReps <- length(fileList[[1]])
  outDepth <- max(length(quantiles),1)


  ##########
  # Initialize Text
  ##########
  if(verbose){
    # what we're doing
    cat("Analyzing", numStage, "directories.\n")

    # whether or not we remove files
    if(rem_file){
      cat("\tRemoving original files.\n\n")
    } else {
      cat("\tNot removing original files.\n\n")
    }
  }


  ##########
  # Dispatch Loop
  ##########
  for(analysis in 1:numStage){

    # indication that something is working
      if(verbose){
        cat("Starting analysis", analysis, "of", numStage, " ...  ")
      }

    # begin logic tree
    if(stage[analysis] %in% c("E","L","P","M","U","FS","FE","FI","F")){
      base_MQ(fList = fileList[[analysis]],oDir = write_dir,sName = stage[analysis],
              nodeNames = mosyNodeNames,nNodes = numMosyNodes,genos = mosyGenos,
              nGenos = numMosyGenos,times = times,nTimes = nTimes,num_repss = numReps,
              mean = mean,quantiles = quantiles,oDepth = outDepth)

    } else if(stage[analysis] == "H"){
      base_MQ(fList = fileList[[analysis]],oDir = write_dir,sName = stage[analysis],
              nodeNames = hNodeNames,nNodes = numHNodes,genos = hGenos,
              nGenos = numHGenos,times = times,nTimes = nTimes,num_repss = numReps,
              mean = mean,quantiles = quantiles,oDepth = outDepth)

    } else if(stage[analysis] %in% c("EE","LE","PE")){
      # get first letter for matching
      holdELP <- nELP[strsplit(x = stage[analysis],split = "",fixed = TRUE)[[1]][1]]

      base_MQ(fList = fileList[[analysis]],oDir = write_dir,sName = stage[analysis],
              nodeNames = mosyNodeNames,nNodes = numMosyNodes,
              genos = 1:holdELP,nGenos = holdELP,
              times = times,nTimes = nTimes,num_repss = numReps,
              mean = mean,quantiles = quantiles,oDepth = outDepth)

    } else if(stage[analysis] == "FSEI"){
      base_MQ(fList = fileList[[analysis]],oDir = write_dir,sName = stage[analysis],
              nodeNames = mosyNodeNames,nNodes = numMosyNodes,
              genos = c("S",paste0("E",1:nEIP),"I"),
              nGenos = nEIP+2,
              times = times,nTimes = nTimes,num_repss = numReps,
              mean = mean,quantiles = quantiles,oDepth = outDepth)
    } # end if tree

    if(verbose){
      cat("Done\n")
    }

  } # end analysis loop

  # check if removing original output
  if(rem_file){file.remove(unlist(fileList))}

} # end function


################################################################################
# Base Statistical Summary Function
################################################################################

#' Base Summary Function
#'
#' This function does the actual calculations for \code{\link{summarize_stats_CSV}}.
#' It calculates mean and quantiles, writing output to the appropriate folder.
#'
#'
#' @param fList File list, all files for this stage, organized by repetition
#' @param oDir Output directory
#' @param sName Stage signifier
#' @param nodeNames Properly formatted vector of node names for printing
#' @param nNodes Number of nodes in the simulation
#' @param genos Vector of genotypes for the header
#' @param nGenos Number of genotypes
#' @param times Vector of sampling times
#' @param nTimes Number of sampled times
#' @param num_repss Number of repetitions from the simulation
#' @param mean Boolean, calculate mean or not
#' @param quantiles Vector of quantiles to calculate, or NULL
#' @param oDepth Max(1, number of quantiles)
#'
#' @return None
#'
#' @importFrom stats quantile
#'
base_MQ <- function(fList,oDir,sName,nodeNames,nNodes,
                    genos,nGenos,times,nTimes,num_repss,
                    mean,quantiles,oDepth){

  # setup input data holder
  inArray <- array(data = 0, dim = c(nTimes, num_repss, nGenos))

  # setup output data holder
  outArray <- array(data = 0, dim = c(nTimes, nGenos+1, oDepth))

  # set time
  outArray[,1,] <- times

  # loop variables
  numCol <- nGenos + 1
  numTot <- numCol*nTimes
  columnNames <- c("Time",genos)

  # loop over all patches and do stats.
  for(node in 1:nNodes){
    # read in all repetitions for this patch
    #  remove time column after reading in.
    for(nR in 1:num_repss){
      inArray[ ,nR, ] <- matrix(data = scan(file = fList[[nR]][node],
                                            what = numeric(), n = numTot,
                                            sep = ",", skip = 1, quiet = TRUE),
                                nrow = nTimes, ncol = numCol, byrow = TRUE)[ ,-1]
    } # end read loop

    # if the user wants the mean, calculate and put out.
    if(mean){
      # loop over genotypes to summarize
      for(whichCol in 1:nGenos){
        outArray[ ,whichCol+1,1] <- .rowMeans(x = inArray[ , ,whichCol], m = nTimes, n = num_repss)
      }

      # write output
      write.table(x = outArray[ , ,1],
                  file = file.path(oDir,file.path(sName,"_Mean_", nodeNames[node], ".csv", fsep = "") ),
                  sep = ",", row.names = FALSE, col.names = columnNames, quote = FALSE)
    } # end mean


    # if the user wants quantiles, do them and write output
    if(!is.null(quantiles)){
      for(whichCol in 1:nGenos){
        outArray[ ,whichCol+1, ] <- t(apply(X = inArray[ , ,whichCol], MARGIN = 1,
                                            FUN = quantile, probs = quantiles,
                                            names = FALSE, type = 8))
      }#end loop to calculate quantiles

      # write output
      for(whichQuant in 1:oDepth){
        # file names
        fName <- file.path(oDir, file.path(sName, "_Quantile_",
                                           formatC(x = quantiles[whichQuant], digits = 4,
                                                   format = "f", decimal.mark = "",
                                                   big.mark = NULL),
                                           "_", nodeNames[node], ".csv", fsep = "") )

        # write output
        write.table(x = outArray[ , ,whichQuant], file = fName, sep = ",",
                    row.names = FALSE, col.names = columnNames, quote = FALSE)
      } # end loop to write files
    } # end quantiles

  } # end loop over nodes

} # end base function
