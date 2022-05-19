################################################################################
#
#   MGDrivE2: Auxiliary data wrangling functions for .csv input/output
#   Marshall Lab
#   Jared B. Bennett (jared_bennett@berkeley.edu)
#   May 2019
#
################################################################################
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

################################################################################
# Plotting and Metrics Summary Function
################################################################################
# This function differs from teh 2 above by being more generic - designed to work
# with V2 and V1 - and by saving data in a ggplot-friendly format.

#' Plotting and Basic Metrics for ggplot2
#'
#' This function is designed to read in all output from one parameter set and
#' then output ggplot-friendly raw data and simple metrics. The raw data maintains
#' the repetition, sex, patch, simulation time, and genotype labels. The metrics are
#' calculated over repetitions and thus maintain sex, patch, simulation time, and
#' genotype labels.
#'
#' This function is generic to \code{MGDrivE} - working on both V1 and V2. Because
#' of this, it has some assumptions about how the data is organized:
#' \itemize{
#'  \item read_dir
#'  \itemize{
#'         \item repetition 1
#'         \itemize{
#'                \item M_[patch1].csv
#'                \item M_[patch2].csv
#'                \item F_[patch1].csv
#'                \item F_[patch2].csv
#'                \item ...
#'         }
#'         \item repetition 2
#'         \itemize{
#'                \item M_[patch1].csv
#'                \item M_[patch2].csv
#'                \item F_[patch1].csv
#'                \item F_[patch2].csv
#'                \item ...
#'         }
#'         \item repetition 3
#'         \item ...
#'  }
#' }
#'
#' This function requires V2 data to be output in V1 format. It will ignore aquatic
#' stages if present. Because of the different way that patches are labeled between
#' versions, it cannot find/check the patches through string manipulation. Therefore,
#' it assumes the shape of the directories and checks that, which indirectly checks
#' the patches, but does not directly check them for correctness or order.
#'
#' Since V1 and V2 use different parameters, and we may be working on data from either,
#' this function reads in an example file to get simulation parameters and check
#' genotypes against.
#'
#' \code{write_dir} defaults to the read directory, placing the analysis with
#' the repeititon data. If this is not desired, the user must build their own
#' \code{write_dir} and supply it to the function. Additionally, the user may
#' supply a name for the analysis using the \code{name} parameter, and all
#' analysis are returned in a folder named \code{analysis_name}.
#'
#' \code{sex} provides the user options on how to handle the male and female
#' data. There are 4 options:
#' \itemize{
#'  \item \code{both}: This reads and keeps male and female data. It is the default
#'  behavior.
#'  \item \code{m} or \code{male}: only analyzes male data. It does not read in female
#'  data or analyze it at all.
#'  \item \code{f} or \code{female}: only analyzes female data. Does not read or
#'  analyze male files.
#'  \item \code{agg}: This reads in both male and female data, but then aggregates
#'  it, so that full population data remains.
#' }
#'
#' \code{goi} can be ignored or provided as a named list. This allows the user to
#' aggregate genotypes, or alleles, in any desired fashion. The default value, \code{NULL},
#' assumes no aggregate or renaming, so all genotypes are used and kept separate.
#' When specific analysis are desired, the input is a named list, e.g.
#' \code{goi=list("A"=c(geno1, geno2), "B"=c(geno1, geno3))}. The list names are
#' returned in the output and can be whatever the user would like to plot. The
#' genotypes within the list must correspond to genotypes in the data and are checked
#' for consistency. This function aggregates the data counts, so if one wishes to
#' do allele frequencies, the genotypes can be duplicated as many times as required
#' (e.g., a diploid would need them duplicated twice, triploid three times, etc.).
#' There is no "total" provided by this function - if a "total" count is desired,
#' then the \code{goi} list must contain \code{goi=list(..."total"=c(all genotypes listed))}.
#'
#' @param read_dir directory with repetition directories in it
#' @param write_dir directory to build output directory, default is \code{read_dir}
#' @param name name for the analysis being performed, a valid string
#' @param sex string dictating how to handle male/female data
#' @param patch_agg boolean, should patches be aggregated, default is \code{FALSE}
#' @param goi named list dictating "genotypes-of-interest", default is NULL
#' @param drop_zero_goi boolean, should genotypes with zero counts be dropped, default is \code{TRUE}
#'
#' @importFrom utils read.csv write.table
#' @importFrom stats quantile
#'
#' @return Builds an analysis directory within \code{write_dir} with three files:
#' \itemize{
#'  \item params.txt: Text file listing the parameters used for the analysis
#'  \item plot.csv: Comma-separated file of simulation data in long-form for ggplot2
#'  \item metrics.csv: Comma-separated file of some basic metrics in long-form for ggplot2
#' }
#'
#' @export
analyze_ggplot_CSV <- function(read_dir, write_dir = read_dir, name = "analysis",
                               sex = "both", patch_agg = FALSE,
                               goi = NULL, drop_zero_goi = TRUE){

  ########################################
  # Parameter Checks
  ########################################
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
  sexString <- paste0('sex: ', paste0(sex, collapse = ', ')) # need the parameter later
  sex <- switch(EXPR = sex, 'm'='male', 'f'='female',
                'both'=c('female','male'), 'agg'=c('female','male'), sex)

  # check patch_agg
  if(!is.logical(patch_agg)){
    stop("'patch_agg' must be 'TRUE' or 'FALSE'")
  }

  # check drop_zero_goi
  if(!is.logical(drop_zero_goi)){
    stop("'drop_zero_goi' must be 'TRUE' or 'FALSE'")
  }


  ########################################
  # Setup
  ########################################
  ##########
  # Files
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
  dirList <- list.dirs(path = read_dir, full.names = TRUE, recursive = FALSE)
  dirList <- grep(pattern = 'analysis_', x = dirList, fixed = TRUE, invert = TRUE, value = TRUE)
  fileList <- lapply(X = dirList,
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

  ##########
  # Variables
  ##########
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
    #  this maintains the shape that GOI needs to be
    goi <- as.list(x = setNames(object = genotypes, nm = genotypes))
  } else if(!all(unlist(x = goi, use.names = FALSE) %in% genotypes)) {
    stop("elements of 'goi' must be valid genotypes in the files, please check goi list")
  }
  numGOI <- length(goi)
  goiNumeric <- lapply(X = goi, FUN = match, table = genotypes)
  goiLengths <- lengths(x = goi)

  ########################################
  # Data
  ########################################
  ####################
  # Setup and Read
  ####################
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


  ####################
  # Organize
  ####################
  ##########
  # Patch Agg.
  ##########
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

  ##########
  # Sex Agg.
  ##########
  if(sexBool){
    # combine over sex
    dataArray <- dataArray[ ,1, , , ,drop = FALSE] + dataArray[ ,2, , , ,drop = FALSE]
    # update dimension name
    dimnames(dataArray)$Sex <- "agg"
  } # end sex reduction

  ##########
  # GOI Reduct.
  ##########
  if(drop_zero_goi){
    # get indices of non-zero GOI
    keepIdx <- which(vapply(X = seq_len(numGOI),
                            FUN = function(x){sum(dataArray[ , , , ,x])>=1},
                            FUN.VALUE = logical(length = 1)))
    # safety check
    if(length(keepIdx)==0){
      stop("There are no genotypes of interest in the analysis.
       Please check simulation output or the elements of 'goi'")
    }
    # subset, handles names internally
    dataArray <- dataArray[ , , , ,keepIdx,drop = FALSE]
  } # end drop check


  ########################################
  # Output and Analysis
  ########################################
  ##########
  # Parameters
  ##########
  aDir <- file.path(write_dir, paste0("analysis_", name))
  dir.create(path = aDir)

  # open connection
  rFileCon <- file(description = file.path(aDir, 'params.txt'), open = 'wt')

  # write parameter names and values
  #  could lapply this and not copy paste?
  writeLines(text = paste0('read_dir: ',read_dir), con = rFileCon, sep = '\n')
  writeLines(text = paste0('write_dir: ',write_dir), con = rFileCon, sep = '\n')
  writeLines(text = sexString, con = rFileCon, sep = '\n')
  writeLines(text = paste0('patch_agg: ', patch_agg), con = rFileCon, sep = '\n')
  # file.path actually combines the list in a really, really nice way!
  #  It does keep extra qutoes, so removed those, then added tabs to make it print nicer
  writeLines(text = paste0('goi: \t', paste0(gsub(pattern = "\"", replacement = "",
                                                  fixed = TRUE, x = file.path(names(goi),
                                                                              goi, fsep = '=')),
                                             collapse = "\n\t")),
             con = rFileCon, sep = '\n')
  writeLines(text = paste0('drop_zero_goi: ', drop_zero_goi), con = rFileCon, sep = '\n')

  # close parameter connection
  close(con = rFileCon)


  ##########
  # Plotting
  ##########
  #  We need to convert wide-form data to long-form (cause ggplot stinks)
  #  if we count well, and replicate labels properly, we can setup the vectors
  #  in one step, without having to add loops or grow objects
  # First, get variables necessary for building the output array.
  # Second, replicate the naming variables
  #  Based on how R does memory, we have to replicate each value equal to the
  #  product of values in previous indexes, then replicate that group by the
  #  product of values in future indexes. This is done with a double-rep() statement
  #   rep( rep(values, each), times)
  #     where "each"is the previous indexes and "times" the future ones
  #     and the first/last indexes have a special case, since you never count yourself
  #     in the product calculation.
  #  setNames is used so that lapply keeps the index names, which data.frame
  #   turns into column names
  # Third, get data. which is super simple (once the naming is done)
  nameDims <- dimnames(dataArray)
  numDims <- length(nameDims)
  numDLen <- lengths(nameDims)

  plotData <- data.frame(
    # state labels
    lapply(X = setNames(object = 1:numDims, nm = names(nameDims)),
           FUN = function(x){
             rep.int(x = rep(x = nameDims[[x]],
                             each = if(x == 1) 1 else prod(numDLen[1:(x-1)]) ),
                     times = if(x == numDims) 1 else prod(numDLen[(x+1):numDims]) )}),
    # value
    "Count" = as.vector(dataArray)
  )

  # write output
  write.table(x = plotData, file = file.path(aDir, 'plot.csv'),
              sep = ',', row.names = FALSE, col.names = TRUE)


  ##########
  # Metrics
  ##########
  # Very similar to above - only difference is that we now collapse over the
  #  "Repetitions" dimension.
  # To simplify this, we calculate the quantile first, then get dimensions from that
  # Lots of the dimensions are hard-coded, and the repetitions must be the first
  #  dimension, in case that changes ever.
  applyDims <- 2:5
  quantData <- apply(X = dataArray, MARGIN = applyDims, FUN = quantile,
                     probs = c(0.025, 0.5, 0.95), type = 8)

  nameDims <- dimnames(quantData)[-1]
  numDims <- length(nameDims)
  numDLen <- lengths(nameDims)

  metricData <- data.frame(
    # set options - didn't like the "%" in the quantiles
    check.names = FALSE,
    # state labels
    lapply(X = setNames(object = 1:numDims, nm = names(nameDims)),
           FUN = function(x){
             rep.int(x = rep(x = nameDims[[x]],
                             each = if(x == 1) 1 else prod(numDLen[1:(x-1)]) ),
                     times = if(x == numDims) 1 else prod(numDLen[(x+1):numDims]) )}),
    # values
    "Min" = as.vector(apply(X = dataArray, MARGIN = applyDims, FUN = min)),
    "Mean" = as.vector(apply(X = dataArray, MARGIN = applyDims[-1], FUN = .colMeans, m=numReps, n=numDLen[['Sex']])),
    "Max" = as.vector(apply(X = dataArray, MARGIN = applyDims, FUN = max)),
    "2.5%" = as.vector(quantData[1, , , , ]),
    "50%" = as.vector(quantData[2, , , , ]),
    "97.5%" = as.vector(quantData[3, , , , ])
  )

  # write output
  write.table(x = metricData, file = file.path(aDir, 'metrics.csv'),
              sep = ',', row.names = FALSE, col.names = TRUE)

} # end analyze_ggplot_CSV

################################################################################
# Plotting and Metrics Summary Function - Helpers
################################################################################
# Taking most of the hashcode stuff from R.oo
#  R has no hashcode ability, and no ability to convert characters to their ASCII
#  values.
#
# Everything below here is right, but doesn't work.
# Numerical issues in the hashcode function returns NAs.
# This is true in R.oo as well, and why they only tested on _short_ strings.
# I'm not sure what a solution is for this.

# ####################
# # ASCII conversion table
# ####################
# # https://github.com/HenrikBengtsson/R.oo/blob/develop/R/ASCII.R
#
# # Idea by Peter Dalgaard, Dept. of Biostatistics, University of Copenhagen, Denmark.
# ASCII <- c("", sapply(1:255, function(i) parse(text=paste("\"\\",
#                    structure(i,class="octmode"), "\"", sep=""))[[1]]) )
#
# # We removed ASCII 0x00, because it represents an empty string in
# # R v2.7.0 (and maybe some earlier version) and in R v2.8.0 we will get
# # a warning.  However, for backward compatibility we will still use it
# # for version prior to R v2.7.0.  See also email from Brian Ripley
# # on 2008-04-23 on this problem.
# if(compareVersion(as.character(getRversion()), "2.7.0") < 0) {
#   ASCII[1] <- eval(parse(text="\"\\000\""))
# }
#
# ####################
# # Int Conversion
# ####################
# # taken from R.oo
# #  not sure this is really the "right" handling - should wrap around and lose
# #  the significant bits.
# # https://en.wikipedia.org/wiki/Integer_overflow
# D2I <- function(x) {
#   intMin <- -2147483648
#   intMax <- 2147483647
#   intRange <- intMax - intMin + 1
#
#   return(as.integer( (x-intMin) %% intRange + intMin ))
# }
#
# ####################
# # Hashcode
# ####################
# # this is loosely based on R.oo "hashcode" function, which is heavily based on
# #  the java string hashcode, which sucks.
# # replacing primes 31 and 0 with 109 and 1
# # vectorizing the calculation
# # https://vanilla-java.github.io/2018/08/12/Why-do-I-think-Stringhash-Code-is-poor.html
# # https://stackoverflow.com/questions/15518418/whats-behind-the-hashcode-method-for-string-in-java
# # https://stackoverflow.com/questions/113511/best-implementation-for-hashcode-method-for-a-collection
# genHash <- function(str2hash){
#
#   # split string into chars
#   strVec <- strsplit(x = str2hash, split = NULL, fixed = TRUE)[[1]]
#
#   # convert chars to ascii values
#   #  this is taken from R.oo
#   intVec <- match(x = strVec, table = ASCII)
#
#   # vectorize hash calculation
#   lenIV <- length(intVec)
#   hashD <- sum(109L^lenIV, 109L^((lenIV-1L):0L) * intVec)
#
#   # convert to integer and return
#   return(D2I(hashD))
# }
