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
#' This function tracks state variables by default; an optional argument \code{Sout}
#' can be provided to track number of event firings each time step (for discrete stochastic simulations),
#' or cumulative intensity (for continuous stochastic simulations), or the rate function of
#' particular events for ODE simulation. The matrix must have number of columns equal to
#' number of events in the system (the number of hazard functions), and a row for each tracking
#' variable. The function \code{\link{track_hinf}} is provided, which builds a matrix to track
#' human infection events.
#'
#' To save output as .csv files, see \code{\link{sim_trajectory_CSV}}.
#'
#'
#' @param x0 the initial marking of the SPN (initial state, M0)
#' @param tmax the final time to end simulation (all simulations start at 0)
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param dt_stoch time-step used for approximation of hazards
#' @param num_reps number of repetitions to run, default is 1.
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param hazards list of hazard functions
#' @param Sout an optional matrix to track event firings
#' @param sampler determines sampling algorithm, one of; "ode", "tau", "cle", or "dm"
#' @param method if \code{sampler} is "ode", the solver to use, from \code{deSolve}
#' @param events a \code{data.frame} of events, may be set to \code{NULL} if not used
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param verbose print a progress bar?
#' @param ... further named arguments passed to the step function
#'
#' @return a list with 2 elements: "state" is the array of returned state values, and "events" will
#'        return events tracked with \code{Sout} if provided, otherwise is \code{NULL}
#'
#' @export
sim_trajectory_R <- function(
  x0, tmax, dt=1, dt_stoch = 0.1, num_reps=1,
  S, hazards, Sout = NULL, sampler = "tau", method = "lsoda",
  events = NULL, batch = NULL, verbose = TRUE, ...){

  if((sampler=="ode") && !is.null(batch)){
    stop("batch migration is incompatible with deterministic simulations")
  }

  # check x0
  base_x0(x0)

  # check/setup step function
  stepFun <- base_stepFunc(
    sampler = sampler,S = S,hazards = hazards,Sout = Sout,
    dt_stoch = dt_stoch,method = method,...
  )

  # check/setup simulation time
  simTimes <- base_time(tt = tmax,dt = dt)

  # check/organize events
  events <- base_events(x0 = x0, events = events, dt = dt)

  # pass everything down to base function
  #  base function does return
  sim_trajectory_base_R(
    x0 = switch(EXPR = sampler, "ode" = x0, round(x0)),
    times = simTimes, num_reps = num_reps, stepFun = stepFun,
    events = events, batch = batch, Sout = Sout, verbose = verbose
  )
}


################################################################################
# Base Functions
################################################################################
#######################################
# Step Funtion
#######################################

# This sets up the step function for sampling
# It's the exact same in sim_trajectory_CSV, so using function
# sample = string denoting type of sampler
# S = stoichiometry matrix
# hazards = list of hazard functions and a flag
# dt_stoch = time-step for approximate hazards
# method = sampler for ODE solver
#
base_stepFunc <- function(sampler,S,hazards,Sout = NULL,dt_stoch,method,...){

  ##########
  # checks for step function!
  ##########
  #  sampler
  if(!(sampler %in% c("ode","tau","cle","dm")) ){
    stop("Sampler must be one of: ode, tau, cle, or dm")
  }

  #  check hazards, depending on sampler
  if((sampler %in% c("ode","cle")) && hazards$flag){
    warning(paste0("For numerical stability, it is strongly recommended that ",
                   "ode and cle samplers use approximate hazards."))
  }

  #  warning to use sparse matrix
  if(attr(class(S),"package") != "Matrix"){
    warning("the stoichiometry 'S' matrix is not sparse.\n significant speed issues possible")
  }


  ##########
  # setup step function!
  #########
  hazFunc <- function(M,t,idx){
    vapply(X = hazards$hazards[idx],FUN = function(h){h(t=t,M=M)}, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
  }

  # generate stepFun
  if(sampler == "tau"){
    stepFun <- step_PTS(S=S, Sout = Sout, haz=hazFunc, sIDX=hazards$approxS, dt = dt_stoch, ...)
  } else if(sampler == "cle"){
    stepFun <- step_CLE(S=S, Sout = Sout, haz=hazFunc, sIDX=hazards$approxS, dt = dt_stoch, ...)
  } else if(sampler == "dm"){
    stepFun <- step_DM(S=S, Sout = Sout, haz=hazFunc, sIDX=hazards$approxS, ...)
  } else if(sampler == "ode"){
    stepFun <- step_ODE(S=S, Sout = Sout, haz=hazFunc, sIDX=hazards$approxS, method = method, ...)
  } else {
    stop("option 'sampler' must be one of 'tau', 'cle', 'dm', 'ode'")
  }

  # return step function
  return(stepFun)
}


#######################################
# Time Function
#######################################

# check time, and setup sampling times
# t0 = initial time to begin simulation
# tt = the final time to end simulation
# dt = the time-step at which to return output
#
base_time <- function(t0 = 0,tt,dt){

  # number of steps we need to take
  n <- (tt-t0) %/% dt + 1
  times <- seq(from=t0,to=tt,by=dt)

  if(length(times) != n){
    stop("error in sequence of times; make sure tt is evenly divisible by dt")
  }

  # return
  return(times)
}


#######################################
# Events Function
#######################################

# place for checks on the initial markings
# doesn't return anything, just errors if there are problems
# x0 = initial petri Net markings
#
base_x0 <- function(x0){
  # check names are valid
  if(any(is.null( names(x0) )) ){
    stop("Marking(s) on x0 is null, please check and try again.")
  }

  # check values are valid
  if(any(x0<0)){
    stop("All states in x0 must be >=0.")
  }

  # any more checks?
}

#######################################
# Events Function
#######################################

# check and organize events
# x0 = the initial marking of the SPN (initial state)
# events = a \code{data.frame} of events
# dt = the time-step at which to return output
#
base_events <- function(x0, events, dt){

  # if no events, escape
  if(is.null(events)) return(NULL)

  # grab state names
  xNames <- names(x0)

  # setup return matrix
  #  hardcoding names so we can reorder them here. May be an issue in upgrades
  retMat <- matrix(data = 0, nrow = nrow(events), ncol = ncol(events),
                   dimnames = list(NULL, c("time","method","var","value")))


  ##########
  # Time
  ##########
  # check event timing, fix if necessary, sort
  # make sure events occur at times that will actually be sampled
  if(any(events$time %% dt != 0)){
    warning("event times do not correspond to a multiple of dt.\n",
            "event times will be rounded up to the nearest time-step!")
    events$time = events$time + (events$time %% dt)
  }

  # get order
  eventOrder <- order(events$time)

  # put times in return matrix
  retMat[ ,"time"] <- events$time[eventOrder]


  ##########
  # Var
  ##########
  # check event name
  #  this is valid for all types of events
  if(!all(events$var %in% xNames)){
    stop("Events ", paste(events$var[!(events$var %in% xNames)],collapse = ", "),
         " are not valid places in the network.\nPlease check the naming.")
  }

  # Convert to numeric indexing,
  # put into return matrix
  retMat[ ,"var"] <- match(x = events$var[eventOrder], table = xNames)


  ##########
  # Method and Value
  ##########
  # check method name
  #  we only support 2 types of events currently
  #  make sure to expand as we add new release types
  if(!all(unique(events$method) %in% c("add","swap"))){
    stop(paste0(unique(events$method)[!(unique(events$method) %in% c("add","swap"))],
                " is not a recognized event method.", collapse = "\n\t"))
  }

  #  we only implemented
  # sort methods and values for method check
  sortMet <- events$method[eventOrder]
  sortVal <- events$value[eventOrder]

  # check release "method"
  #  turn method strings into numeric values
  #  check value, as it depends on method

  # "add" method
  addIdx <- which(sortMet == "add")
  if(length(addIdx) > 0){
    # convert method to numeric
    # add to return matrix
    retMat[addIdx,"method"] <- 1

    # put values into return matrix
    retMat[addIdx,"value"] <- as.numeric(sortVal[addIdx])
  }

  # swap method
  swapIdx <- which(sortMet == "swap")
  if(length(swapIdx) > 0){
    # convert method to numeric
    # add to return matrix
    retMat[swapIdx,"method"] <- 2

    # check values for proper naming
    if(!all(sortVal[swapIdx] %in% xNames)){
      stop("Swap events ", paste(sortVal[swapIdx][!(sortVal[swapIdx] %in% xNames)],collapse = ", "),
           " are not valid places in the network.\nPlease check the values.")
    }

    # Convert to numeric indexing,
    # put into return matrix
    retMat[swapIdx,"value"] <- match(x = sortVal[swapIdx], table = xNames)
  }

  # next method addition
  #  expand release types here, follow same patters


  # return updated events
  return(retMat)
}


################################################################################
# Base simulation trajectory
################################################################################

#' Simulate Trajectory From one SPN Model
#'
#' This is an internal function to \code{\link{sim_trajectory_R}}. It does the
#' actual sampling once all of the functions have been checked and setup.
#'
#'
#' @param x0 the initial marking of the SPN (initial state)
#' @param times sequence of sampling times
#' @param num_reps number of repetitions to run
#' @param stepFun a sampling function
#' @param events a \code{data.frame} of events (uses the same format as required in package \code{deSolve} for consistency, see \code{\link[deSolve]{events}} for more information)
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param Sout an optional matrix to track event firings
#' @param verbose print a progress bar?
#'
#' @return matrix of sampled values
#'
#' @importFrom stats rbinom
sim_trajectory_base_R <- function(x0, times, num_reps, stepFun, events = NULL,
                                  batch = NULL, Sout = NULL, verbose = TRUE){

  # setup return array
  nTime <- length(times)
  retArray <- array(data = 0, dim = c(nTime, length(x0)+1, num_reps),
                    dimnames = list(times, c("time",names(x0)), 1:num_reps))
  retArray[ ,1, ] <- times
  retArray[1,-1, ] <- x0

  # set up event tracking (tracks differences, so one less time output than state)
  if(!is.null(Sout)){
    track <- TRUE
    ret_events <- array(data = 0, dim = c(nTime-1, nrow(Sout)+1, num_reps),
                        dimnames = list(times[-1], c("time",rownames(Sout)), 1:num_reps) )
    ret_events[,1,] <- times[-1]
  } else {
    track <- FALSE
    ret_events <- NULL
  }

  # setup event info
  eventsNotNull <- !is.null(events)
  eventIdx <- 1L
  eventLen <- nrow(events)

  # setup batch migration info
  batchNotNull <- !is.null(batch)
  batchIdx <- 1L
  batchLen <- length(batch)


  # loop over num_reps, calling base function
  for(r in 1:num_reps){

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
      retArray[i,-1,r] <- state$x
      if(track){
        ret_events[i-1,-1,r] <- state$o
      }

      # progress bar
      if(verbose){setTxtProgressBar(pb = pbar,value = i)}
    } # end sim loop

    if(verbose){
      close(pbar)
      cat(" --- end simulation --- \n")
    }

  } # end repetition loop

  # return stuff
  return(list("state"=retArray,"events"=ret_events))
}
