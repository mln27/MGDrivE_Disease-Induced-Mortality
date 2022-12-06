################################################################################
#
#   MGDrivE2: decoupled sampling trajectory interface
#   Marshall Lab
#   Agastya Mondal (agastya_mondal@berkeley.edu )
#   October 2021
#
################################################################################

################################################################################
# User-facing simulation trajectory
################################################################################

#' Simulate Trajectory From a Decoupled SPN Model
#'
#' This function provides a unified interface to the various decoupled simulation algorithms
#' for SPN, returning output sampled at a lattice of time points to the user, and
#' handling various exogenous events that may occur during the simulation
#' (such as release of adult mosquitoes).
#'
#' \code{dt_stoch} is used by the Poisson Time-Step (\code{\link{step_PTS_decoupled}})
#' method to approximate the hazards. A smaller \code{dt_stoch} provides a better
#' approximation, but will take longer to run.
#'
#' The stoichiometry matrix (\code{S}) is generated in \code{\link{spn_S}}.
#'
#' The list of hazards (\code{hazards}) come from \code{\link{spn_hazards}}.
#'
#' Two samplers are provided. The default is a Poisson Time-Step
#' (\code{\link{step_PTS_decoupled}}) method.
#' Additionally, for convenience, an ODE "sampler" (\code{\link{step_ODE_decoupled}}) is
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
#' To save output as .csv files, see \code{\link{sim_trajectory_CSV_decoupled}}.
#'
#'
#' @param x0 the initial marking of the SPN (initial state, M0)
#' @param h0 the initial human state distribution
#' @param inf_labels labels corresponding to female mosquito infection hazard
#' @param tmax the final time to end simulation (all simulations start at 0)
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param dt_stoch time-step used for approximation of hazards
#' @param num_reps number of repetitions to run, default is 1.
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param hazards list of hazard functions
#' @param Sout an optional matrix to track event firings
#' @param sampler determines sampling algorithm, one of; "tau-decoupled", "ode-decoupled"
#' @param method if \code{sampler} is "ode-decoupled", the solver to use, from \code{deSolve}
#' @param events a \code{data.frame} of events, may be set to \code{NULL} if not used
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param verbose print a progress bar?
#' @param SPN_P stochastic Petri Net places
#' @param theta parameters
#' @param human_ode human ode function
#' @param cube inheritance cube
#' @param h0 the initial human state distribution
#' @param ... further named arguments passed to the step function
#'
#' @return a list with 2 elements: "state" is the array of returned state values, and "events" will
#'        return events tracked with \code{Sout} if provided, otherwise is \code{NULL}
#'
#' @export
sim_trajectory_R_decoupled <- function(
    x0, h0, tmax, inf_labels, dt = 1, dt_stoch = 0.1, num_reps = 1,
    S, hazards, SPN_P, theta, Sout = NULL, sampler = "tau", method = "lsoda",
    events = NULL, batch = NULL, verbose = TRUE, human_ode = "Imperial",
    cube, ...){

  if(sampler == "ode-decoupled" & !is.null(batch)){
    stop("batch migration is incompatible with deterministic simulations")
  }

  # check x0
  base_x0(x0)

  # check/setup step function
  stepFun <- base_stepFunc_decoupled(
    sampler = sampler, S = S, hazards = hazards, Sout = Sout,
    dt_stoch = dt_stoch, method = method, inf_labels = inf_labels,
    human_ode = human_ode, ... )

  # check/setup simulation time
  simTimes <- base_time(tt = tmax, dt = dt)

  # check/organize events and x0
  events <- base_events(x0 = x0, events = events, dt = dt)

  if(human_ode == "SIS"){
    sampling_func <- sim_trajectory_base_R_decoupled_SIS
  } else if(human_ode == "Imperial"){
    sampling_func <- sim_trajectory_base_R_decoupled_Imperial
  } else {
    stop("Only SIS or Imperial models are currently supported for human_ode")
  }

  # pass everything down to base function
  #  base function does return
  sampling_func(
    x0 = switch(EXPR = sampler, "ode-decoupled" = x0, round(x0)),
    h0 = h0, SPN_P = SPN_P, theta = theta, times = simTimes,
    num_reps = num_reps, stepFun = stepFun, events = events,
    batch = batch, Sout = Sout, verbose = verbose, cube = cube
  )
}

################################################################################
# Base Functions
################################################################################
#######################################
# Step Function
#######################################

# This sets up the step function for sampling
# It's the exact same in sim_trajectory_CSV, so using function
# sample = string denoting type of sampler
# S = stoichiometry matrix
# hazards = list of hazard functions and a flag
# dt_stoch = time-step for approximate hazards
# method = sampler for ODE solver
# inf_labels = ??
# human_ode = ??
#
base_stepFunc_decoupled <- function(sampler,S,hazards,Sout = NULL, dt_stoch,
                                    method,inf_labels,human_ode,...){
  ##########
  # checks for step function!
  ##########
  #  sampler
  # for decoupled sampling we only have two samplers for now
  # tau (stochastic mosquito sampler) and ode (deterministic mosquito sampler)
  if(!(sampler %in% c("tau-decoupled", "ode-decoupled"))){
    stop("Sampler must be one of: ode-decoupled or tau-decoupled. Other samplers not yet supported.")
  }

  #  check hazards, depending on sampler
  if((sampler %in% c("ode-decoupled")) && hazards$flag){
    warning(paste0("For numerical stability, it is strongly recommended that ",
                   "ode and cle samplers use approximate hazards."))
  }

  #  warning to use sparse matrix
  if(attr(class(S), "package") != "Matrix"){
    warning("the stoichiometry 'S' matrix is not sparse.\n significant speed issues possible")
  }


  ##########
  # setup step function! Imperial requires past values via human_trace
  #########
  if(human_ode == "Imperial") {

    hazFunc <- function(M, t, h, human_trace = NULL, idx) {
      vapply(X = names(hazards$hazards[idx]), FUN = function(transition_name){
        # we need to pass the human state and past values (trace) to only the female infection
        # hazards, but not the rest
        haz_closure <- hazards$hazards[[transition_name]]
        if(transition_name %in% inf_labels){
          return(haz_closure(t=t,M=M,h=h,human_trace=human_trace))
        } else {
          return(haz_closure(t=t,M=M))
        }
      }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
    }

  } else {

    hazFunc <- function(M, t, h, human_trace = NULL, idx) {
      vapply(X = names(hazards$hazards[idx]), FUN = function(transition_name){
        # we need to pass the human state to only the female infection
        # hazards, but not the rest
        haz_closure <- hazards$hazards[[transition_name]]
        if (transition_name %in% inf_labels) {
          return(haz_closure(t = t, M = M, h = h))
        } else {
          return(haz_closure(t = t, M = M))
        }
      }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
    }

  }


  # generate stepFun
  if(sampler == "tau-decoupled"){
    stepFun <- step_PTS_decoupled(S=S, Sout=Sout, haz=hazFunc, sIDX=hazards$approxS,
                                  dt=dt_stoch, human_ode=human_ode, ...)
  } else if(sampler == "ode-decoupled"){
    stepFun <-step_ODE_decoupled(S=S, Sout=Sout, haz=hazFunc, sIDX=hazards$approxS,
                                 method=method, human_ode=human_ode, ...)
  } else {
    stop("option 'sampler' must be one of 'tau-decoupled, or 'ode-decoupled'")
  }

  # return step function
  return(stepFun)
}

################################################################################
# Base simulation trajectory
################################################################################

#' Simulate Trajectory From one SPN Model using Imperial Malaria Model
#'
#' This is an internal function to \code{\link{sim_trajectory_R_decoupled}}. It does the
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
#' @param h0 initial human state distribution
#' @param SPN_P stochastic Petri Net, places
#' @param theta parameters
#' @param cube inheritance cube
#'
#' @return matrix of sampled values
#'
#' @importFrom deSolve ode
#' @importFrom stats rbinom
sim_trajectory_base_R_decoupled_Imperial <- function(
    x0, h0, SPN_P, theta, times, num_reps, stepFun,
    events = NULL, batch = NULL, Sout = NULL, verbose = TRUE, cube = NULL){

  # setup return array
  nTime <- length(times)
  na <- theta$na #number of age compartments
  clin_inc <- theta$clin_inc
  human_state_labels <- generate_Imperial_human_state_labels(na)
  retArray <- array(data = 0, dim = c(nTime, length(c(x0, h0, clin_inc)) + 1, num_reps),
                    dimnames = list(times, c("time", c(names(x0), human_state_labels)), 1:num_reps) )
  retArray[, 1, ] <- times
  retArray[1,-1, ] <- c(x0, h0, clin_inc)

  # set up event tracking (tracks differences, so one less time output than state)
  if(!is.null(Sout)){
    track <- TRUE
    ret_events <- array(data = 0, dim = c(nTime - 1, nrow(Sout) + 1, num_reps),
                        dimnames = list(times[-1], c("time", rownames(Sout)), 1:num_reps) )
    ret_events[, 1,] <- times[-1]
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

    # initialize human trace list
    names <- as.character(seq(1, 15))
    human_trace <- vector("list", length(names))
    names(human_trace) <- names

    human_trace[["1"]] <- h0

    state <- list("x" = NULL, "o" = NULL, "h" = NULL)
    state$x <- x0
    state$h <- h0

    eventIdx <- 1L
    batchIdx <- 1L


    # main simulation loop
    for(i in 2:nTime){

      # iterate the step function until this delta t is over
      t0 <- times[i - 1]
      t1 <- times[i]
      dt <- t1 - t0

      state <- stepFun(state$x, state$h, t0, dt, human_trace)

      mosquito_counts <- aggregate_female_SEI(SPN_P, state$x)

      func <- human_Imperial_ODE
      # attach distribution of infectious mosquitoes across genotypes
      # and total number of adult mosquitoes
      theta$I_V <- mosquito_counts$I_V
      theta$N_V <- sum(mosquito_counts$N_V)

      out <- ode(y=state$h, times=c(t0,t1), func=func, parms=theta, rtol=1e-10)

      # update human state matrix
      y <- tail(out, 1)
      num_states <- 10
      # keep clinical incidence separate, it doesn't contribute to ODE dynamics
      # we're just using it for presentation/viz later on
      human_state_matrix <- matrix(y[2:length(y)], ncol = num_states)
      clin_inc_idx <- 10
      clin_inc <- human_state_matrix[, clin_inc_idx]
      human_state_matrix <- human_state_matrix[,-clin_inc_idx]
      colnames(human_state_matrix) <- c("S", "T", "D", "A", "U", "P", "ICA", "IB", "ID")

      state$h <- human_state_matrix
      idx <- as.character(t1 + 1)
      human_trace[[idx]] <- human_state_matrix

      # only keep the last 15 values in memory, otherwise program becomes too slow
      max_entries <- 15
      if(length(human_trace) > max_entries){
        num_entries <- length(human_trace)
        human_trace <- human_trace[(num_entries - 15):num_entries]
      }

      if( all(fequal(state$x, 0)) ){
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

      # record output + add colnames
      retArray[i,-1, r] <- c(state$x, state$h, clin_inc)
      if(track){
        ret_events[i - 1,-1, r] <- state$o
      }

      # progress bar
      if(verbose){setTxtProgressBar(pb = pbar, value = i)}
    } # end sim loop

    if(verbose){
      close(pbar)
      cat(" --- end simulation --- \n")
    }

  } # end repetition loop

  # return stuff
  return(list("state" = retArray, "events" = ret_events))
}


#' Simulate Trajectory From one SPN Model using Human SIS model
#'
#' This is an internal function to \code{\link{sim_trajectory_R_decoupled}}. It does the
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
#' @param h0 initial human state distribution
#' @param SPN_P stochastic Petri Net, places
#' @param theta parameters
#' @param cube inheritance cube
#'
#' @return matrix of sampled values
#'
#' @importFrom deSolve ode
sim_trajectory_base_R_decoupled_SIS <- function(
    x0, h0, SPN_P, theta, times, num_reps, stepFun,
    events = NULL, batch = NULL, Sout = NULL, verbose = TRUE, cube = NULL){

  # setup return array
  nTime <- length(times)
  retArray <- array(data = 0, dim = c(nTime, length(c(x0, h0)) + 1, num_reps),
                    dimnames = list(times, c("time", names(c(x0, h0))), 1:num_reps) )
  retArray[, 1, ] <- times
  retArray[1,-1, ] <- c(x0, h0)

  # set up event tracking (tracks differences, so one less time output than state)
  if(!is.null(Sout)){
    track <- TRUE
    ret_events <- array(data = 0,dim = c(nTime - 1, nrow(Sout) + 1, num_reps),
                        dimnames = list(times[-1], c("time", rownames(Sout)), 1:num_reps) )
    ret_events[, 1,] <- times[-1]
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
      pbar <- txtProgressBar(min = 1, max = nTime, style = 3)
      cat(" --- begin simulation --- \n")
    }

    # reset at the beginning of every simulation
    state <- list("x" = NULL, "o" = NULL, "h" = NULL)
    state$x <- x0
    state$h <- h0

    eventIdx <- 1L
    batchIdx <- 1L


    # main simulation loop
    for(i in 2:nTime){

      # iterate the step function until this delta t is over
      t0 <- times[i - 1]
      t1 <- times[i]
      dt <- t1 - t0

      state <- stepFun(state$x, state$h, t0, dt)

      mosquito_counts <- aggregate_female_SEI(SPN_P, state$x)
      func <- human_SIS_ODE

      out <- as.data.frame(ode(y = state$h, times = c(t0, t1), func = func,
                               parms = list(I_V = mosquito_counts$I_V,
                                            N_V = mosquito_counts$N_V,
                                            a = theta$a, b = cube$b, r = theta$r)
      ))

      # update human state
      state$h <- c(H_S = tail(out$H_S, 1), H_I = tail(out$H_I, 1))

      if( all(fequal(state$x,0)) ){
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
      retArray[i,-1, r] <- c(state$x, state$h)
      if(track){
        ret_events[i - 1,-1, r] <- state$o
      }

      # progress bar
      if(verbose){setTxtProgressBar(pb = pbar, value = i)}
    } # end sim loop

    if(verbose){
      close(pbar)
      cat(" --- end simulation --- \n")
    }

  } # end repetition loop

  # return stuff
  return(list("state" = retArray, "events" = ret_events))
}
