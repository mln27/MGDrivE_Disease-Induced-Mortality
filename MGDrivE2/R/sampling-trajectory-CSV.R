################################################################################
#
#   MGDrivE2: trajectory interface
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# User-facing simulation trajectory
################################################################################

#' Simulate Trajectory From a SPN Model
#'
#' This function provides a unified interface to the various simulation algorithms
#' for SPN, returning output sampled at a lattice of time points to the user, and
#' handling various exogenous events that may occur during the simulation
#' (such as release of adult mosquitoes).
#'
#' \code{dt_stoch} is used by the Poisson Time-Step (\code{\link{step_PTS}}) and
#' Chemical Langevin (\code{\link{step_CLE}}) methods to approximate the hazards.
#' A smaller \code{dt_stoch} provides a better approximation, but will take longer
#' to run.
#'
#' The stoichiometry matrix (\code{S}) is generated in \code{\link{spn_S}}.
#'
#' The list of hazards (\code{hazards}) come from \code{\link{spn_hazards}}.
#'
#' Several samplers are provided. The default is a Poisson Time-Step
#' (\code{\link{step_PTS}}) method. Other options are Gillespie's Direct Method
#' (\code{\link{step_DM}}) and a Chemical Langevin sampler (\code{\link{step_CLE}}).
#' Additionally, for convenience, an ODE "sampler" (\code{\link{step_ODE}}) is
#' provided for compatibility with other samplers. This function uses methods from
#' \code{deSolve}.
#'
#' If using the \code{ode} sampler, several \code{methods} are provided in the \code{deSolve}
#' package, see \code{\link[deSolve]{ode}}. For inhomogeneous systems, consider
#' using the "rk4" method to avoid excessive integration times.
#'
#' Additionally, \code{events} objects may be used for "releases", and must follow
#' the format required by \code{deSolve}. This means that events must be a \code{data.frame}
#' with four columns - \code{var}, \code{time}, \code{value}, and \code{method} -
#' in any order and \code{stringsAsFactors = FALSE}. \code{Events} currently must be one of two types:
#' \itemize{
#'  \item: 'add': When the \code{method} is \code{add}, this represents a classic
#'         release scenario. Appropriate values from the \code{value} column are
#'         added to the state defined in the \code{var} column at the appropriate
#'         moments as defined in the \code{time} column.
#'  \item: 'swap': The \code{swap} method is similar to a small-molecule spray, where
#'         a chemical is applied to the landscape and "turns-on" or "turns-off"
#'         genotypes. This functionality uses the states defined in the \code{value}
#'         column to move the counts within those states to states defined in the
#'         \code{var} column. Then, the state specified in the \code{value} column
#'         is set to \code{0}, effectively "swapping" \code{value} into \code{var}.
#'         Again, these events occur as specified in the \code{time} column.
#' }
#'
#' This function writes all output to .csv files. Each simulation is written to
#' a \code{folder} element - the number of repetitions is the number of folders
#' provided. What life-stages get recorded is specified by the \code{stage} parameter.
#' All life-stages can be stored, or any subset thereof. Females are split by
#' infection status, i.e. by "S", "E", or "I".
#'
#' This function tracks state variables specified by argument \code{stage} by default; an optional argument \code{Sout}
#' can be provided to track number of event firings each time step (for discrete stochastic simulations),
#' or cumulative intensity (for continuous stochastic simulations), or the rate function of
#' particular events for ODE simulation. The matrix must have number of columns equal to
#' number of events in the system (the number of hazard functions), and a row for each tracking
#' variable. If \code{Sout} is provided, it output an additional csv, "events.csv".
#' The function \code{\link{track_hinf}} is provided, which builds a matrix to track
#' human infection events.
#'
#' To return simulations to R for further processing, see \code{\link{sim_trajectory_R}}.
#'
#'
#' @param x0 the initial marking of the SPN (initial state, M0)
#' @param tmax the final time to end simulation
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param dt_stoch time-step used for approximation of hazards
#' @param folders vector of folders to write output
#' @param stage life-stages to print. Any combination of: "E", "L", "P"," M", "U", "F", "H"
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param hazards list of hazard functions
#' @param Sout an optional matrix to track event firings
#' @param sampler determines sampling algorithm, one of; "ode", "tau", "cle", or "dm"
#' @param method if \code{sampler} is "ode", the solver to use, from \code{deSolve}
#' @param events a \code{data.frame} of events
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param verbose print a progress bar?
#' @param ... further named arguments passed to the step function
#'
#' @return NULL - prints output to .csv files
#'
#' @importFrom stats rbinom
#'
#' @export
sim_trajectory_CSV <- function(
  x0, tmax, dt=1, dt_stoch = 0.1, folders="./",
  stage=c("M","F"), S, hazards, Sout = NULL, sampler = "tau",
  method = "lsoda", events = NULL, batch = NULL, verbose = TRUE, ...
){

  if(sampler == "ode" & !is.null(batch)){
    stop("batch migration is incompatible with deterministic simulations")
  }

  # check x0
  base_x0(x0)

  # check/setup step function
  stepFun <- base_stepFunc(
    sampler = sampler,S = S,hazards = hazards, Sout = Sout,
    dt_stoch = dt_stoch,method = method,...
  )

  # check/setup simulation time
  simTimes <- base_time(tt = tmax,dt = dt)

  # check/organize events
  events <- base_events(x0 = x0, events = events, dt = dt)

  # ensure printing only states that exist
  if(!all(stage %in% c("E","L","P","M","U","F","H"))){
    stop("Print stages must be one or more of 'E', 'L', 'P', 'M', 'U', 'F', 'H'.")
  }
  # remove duplicates
  stage <- unique(stage)


  # pass everything down to base function
  sim_trajectory_base_CSV(
    x0 = switch(EXPR = sampler, "ode" = x0, round(x0)),
    times = simTimes, stepFun = stepFun,
    folders = folders, stage = stage,
    events = events, batch = batch, Sout = Sout, verbose = verbose
  )

  # no return
}


################################################################################
# Base simulation trajectory
################################################################################

#' Simulate Trajectory From one  SPN Model
#'
#' This is an internal function to \code{\link{sim_trajectory_CSV}}. It does the
#' actual sampling once all of the functions have been checked and setup.
#'
#' @param x0 the initial marking of the SPN (initial state)
#' @param times sequence of sampling times
#' @param stepFun a sampling function
#' @param folders vector of folders to write output
#' @param stage vector of life-stages to print
#' @param events a \code{data.frame} of events (uses the same format as required
#' in package \code{deSolve} for consistency, see \code{\link[deSolve]{events}}
#' for more information)
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param Sout an optional matrix to track event firings
#' @param verbose print a progress bar?
#'
#' @return no return, prints .csv files into provided folders
#'
#' @importFrom stats rbinom
sim_trajectory_base_CSV <- function(
  x0, times, stepFun, folders,
  stage, events=NULL, batch = NULL,
  Sout = NULL, verbose=TRUE
){

  # stuff for later
  xNames <- names(x0)
  nTime <- length(times)

  # make index list for all possible output
  #  function defined below
  pStuff <- base_print_stuff(stage = stage, state_names = xNames)
  pList <- pStuff$pList
  stage <- pStuff$stage
  lenPS <- length(stage)

  # track rates/events?
  if(!is.null(Sout)){
    track <- TRUE
  } else {
    track <- FALSE
  }

  # setup event info
  eventsNotNull <- !is.null(events)
  eventIdx <- 1L
  eventLen <- nrow(events)

  # setup batch migration info
  batchNotNull <- !is.null(batch)
  batchIdx <- 1L
  batchLen <- length(batch)


  # loop over repetitions
  for(num_reps in 1:length(folders)){

    # progress bar
    if(verbose){
      pbar <- txtProgressBar(min = 1,max = nTime,style = 3)
      cat(" --- begin simulation --- \n")
    }

    # reset at the beginning of every simulation
    state <- list("x"=NULL,"o"=NULL)
    state$x <- x0
    eventIdx <- 1L
    batchIdx <- 1L

    # Initialize files
    fileCons <- setNames(object = vector(mode = "list", length = lenPS), nm = stage)
    for(curS in 1:lenPS){
      # open file connection
      fileCons[[curS]] <- file(description = paste0(folders[num_reps],.Platform$file.sep,stage[curS],".csv"),
                               open = "wt")
      # write header
      writeLines(text = paste0(c("Time",xNames[pList[[curS]] ]),collapse = ","),
                 con = fileCons[[curS]], sep = "\n")
      # write T0
      writeLines(text = paste0(c(0,x0[pList[[curS]] ]),collapse = ","),
                 con = fileCons[[curS]], sep = "\n")
    }

    # tracking rates/event firing
    if(track){
      event_con <- file(
        description = paste0(folders[num_reps],.Platform$file.sep,"Events.csv"),
        open = "wt"
      )
      writeLines(
        text = paste0(c("Time",rownames(Sout)),collapse = ","),
        con = event_con, sep = "\n"
      )
    }


    # main simulation loop
    for(i in 2:nTime){

      # iterate the step function until this delta t is over
      t0 <- times[i-1]
      t1 <- times[i]
      dt <- t1-t0
      # tNow <- times[i]
      state <- stepFun(state$x,t0,dt)

      if(all(fequal(state$x,0))){
        if(verbose){close(pbar)}

        # fill rest of file with 0s
        #  this makes sure all files are the same size/shape for analysis
        for(nT in i:nTime){
          for(curS in 1:lenPS){
            writeLines(text = paste0(c(times[nT],state$x[pList[[curS]] ]),collapse = ","),
                       con = fileCons[[curS]], sep = "\n")
          }
          # tracking rates/event firing
          if(track){
            writeLines(
              text = paste0(c(times[nT],rep(0,nrow(Sout))),collapse = ","),
              con = event_con, sep = "\n"
            )
          }
        } # end print loop

        warning(" --- marking of net is zero; terminating simulation early --- \n")
        break
      }

      # add the event to the state vector
      #  check that there are events
      if(eventsNotNull){
        # check that we haven't performed all the events
        # if there's more events to perform, see if it's the correct time
        while((eventIdx <= eventLen) && (events[eventIdx,"time"] <= t1)){
          # switch on release type
          #  can't name them, cause R is weird. It won't accept numeric labels
          #  also, numeric indices can't have a default, but at least R doesn't cascade
          switch(events[eventIdx,"method"],
                 {# add
                   state$x[events[eventIdx,"var"]] <- state$x[events[eventIdx,"var"]] + events[eventIdx,"value"]
                 },
                 {# swap
                   state$x[events[eventIdx,"var"]] <- state$x[events[eventIdx,"var"]] + state$x[events[eventIdx,"value"]]
                   state$x[events[eventIdx,"value"]] <- 0
                 }
          ) # end switch

          # increment index
          eventIdx <- eventIdx + 1L
        } # end while
      } # end event null check

      # perform batch migration
      #  check that there is batch migration
      if(batchNotNull){
        # check that we haven't performed all the events
        # if there's more events to perform, see if it's the correct time
        while((batchIdx <= batchLen) && (batch[[batchIdx]]$time <= t1)){
          # how many go?
          moving <- rbinom(n = length(batch[[batchIdx]]$from),
                           size = as.integer(state$x[batch[[batchIdx]]$from]),
                           prob = batch[[batchIdx]]$prob)
          # update the state
          state$x[batch[[batchIdx]]$from] <- state$x[batch[[batchIdx]]$from] - moving
          state$x[batch[[batchIdx]]$to] <- state$x[batch[[batchIdx]]$to] + moving

          # increment index
          batchIdx <- batchIdx + 1L
        } # end while
      } # end batch null check

      # record output
      for(curS in 1:lenPS){
        writeLines(text = paste0(c(t1,state$x[pList[[curS]] ]),collapse = ","),
                   con = fileCons[[curS]], sep = "\n")
      }
      # tracking rates/event firing
      if(track){
        writeLines(
          text = paste0(c(t1,state$o),collapse = ","),
          con = event_con, sep = "\n"
        )
      }

      # progress bar
      if(verbose){ setTxtProgressBar(pb = pbar,value = i) }
    } # end sim loop

    # close connections
    for(curS in 1:lenPS){ close(con = fileCons[[curS]]) }
    if(track){
      close(con=event_con)
    }

    if(verbose){
      close(pbar)
      cat(" --- end simulation --- \n")
    }
  } # end repetitions loop

  # no return
}


################################################################################
# Base print
################################################################################

# there was a lot to make sure that everything gets printed properly,
#  so moving it down here
base_print_stuff <- function(stage,state_names){
  # make index list for all possible output
  pList <- lapply(X = stage, FUN = function(x){grep(pattern = paste0("^",x),
                                                           x = state_names)})

  # check for silly user request
  nullIDX <- lengths(pList) == 0
  # warn them about it
  if(any(nullIDX)){
    warning(paste0(stage[nullIDX], " is not present in the sim, removing from output request.\n"))
  }

  # remove it
  stage <- stage[!nullIDX]
  pList[nullIDX] <- NULL

  # check if females are request, split into epi stages if so
  if("F" %in% stage){
    # get index of female
    fIDX <- match(x = "F", table = stage)

    # list of split female names
    fSplitList <- strsplit(x = state_names[pList[[fIDX]] ],
                           split = "_", fixed = TRUE)

    # mosquito-only case
    if(is.na(fSplitList[[1]][4]) || (fSplitList[[1]][4] != "S")){
      # just change female name to match SEI output later
      stage[stage == "F"] <- "FS"

    } else {
      # get index of female, remove from stage
      #  hold on to pList so I can use it
      stage <- stage[-fIDX]


      # get list of female SEI indices
      #  unlist states, then grep them (can't match because E states)
      #  then return index as list
      infStat <- vapply(X = fSplitList, FUN = '[[', 4, FUN.VALUE = character(length = 1))
      fList <- lapply(X = c("S", "E", "I"), FUN = function(x){
        pList[[fIDX]][grep(pattern = x, x = infStat, fixed = TRUE)]
      })

      # now remove old female indices from pList
      pList[fIDX] <- NULL

      # append things and return updated pStages and pList
      stage <- c(stage, "FS", "FE", "FI")
      pList <- c(pList,fList)

    } # end update for F_S/E/I
  } # end check if F

  # return 2 things in list
  return(list("stage" = stage,
              "pList" = pList))
}
