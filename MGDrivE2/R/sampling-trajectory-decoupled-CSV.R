################################################################################
#
#   MGDrivE2: trajectory interface for decoupled sampling
#   Marshall Lab
#   Agastya Mondal (agastya_mondal@berkeley.edu)
#   September 2021
#
################################################################################

################################################################################
# User-facing simulation trajectory
################################################################################

#' Simulate Trajectory From a Decoupled SPN Using the Imperial Malaria Model
#'
#' This function provides a unified interface to the various decoupled simulation algorithms
#' for SPN, returning output sampled at a lattice of time points to the user, and
#' handling various exogenous events that may occur during the simulation
#' (such as release of adult mosquitoes). This function is used in decoupled sampling, where
#' the mosquito and human states are separated. This is used \strong{only} when using the
#' Imperial model of malaria transmission.
#'
#' \code{dt_stoch} is used by the Poisson Time-Step (\code{\link{step_PTS_decoupled}})
#' method to approximate the hazards.
#' A smaller \code{dt_stoch} provides a better approximation, but will take longer
#' to run.
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
#' This function writes all output to .csv files. Each simulation is written to
#' a \code{folder} element - the number of repetitions is the number of folders
#' provided. For now, only adult mosquito states, human states, clinical incidence, and pathogen prevalence are
#' written to CSVs.
#'
#' This function tracks state variables specified by argument \code{stage} by default; an optional argument \code{Sout}
#' can be provided to track number of event firings each time step (for discrete stochastic simulations),
#' or cumulative intensity (for continuous stochastic simulations), or the rate function of
#' particular events for ODE simulation. The matrix must have number of columns equal to
#' number of events in the system (the number of hazard functions), and a row for each tracking
#' variable. If \code{Sout} is provided, it outputs an additional csv, "Tracking.csv".
#' The function \code{\link{track_hinf}} is provided, which builds a matrix to track
#' human infection events.
#'
#' To return simulations to R for further processing, see \code{\link{sim_trajectory_R_decoupled}}.
#'
#'
#' @param x0 the initial marking of the SPN (initial state, M0)
#' @param h0 the initial human state distribution
#' @param inf_labels labels corresponding to female mosquito infection hazard
#' @param tmax the final time to end simulation
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param dt_stoch time-step used for approximation of hazards
#' @param folders vector of folders to write output
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param hazards list of hazard functions
#' @param Sout an optional matrix to track event firings
#' @param sampler determines sampling algorithm, one of; "tau-decoupled", "ode-decoupled"
#' @param method if \code{sampler} is "ode-decoupled", the solver to use, from \code{deSolve}
#' @param events a \code{data.frame} of events
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param verbose print a progress bar?
#' @param SPN_P stochastic Petri Net places
#' @param theta parameters
#' @param human_ode human ode function
#' @param cube inheritance cube
#' @param ... further named arguments passed to the step function
#'
#' @return NULL - prints output to .csv files
#'
#' @export
sim_trajectory_CSV_decoupled <- function(
    x0, h0, inf_labels, tmax, dt = 1, dt_stoch = 0.1, folders = "./",
    S, hazards, SPN_P, theta, Sout = NULL, sampler = "tau", method = "lsoda",
    events = NULL, batch = NULL, verbose = TRUE, human_ode = "Imperial",
    cube, ...){

  if((sampler=="ode-decoupled") && !is.null(batch)){
    stop("batch migration is incompatible with deterministic simulations")
  }

  if(human_ode != "Imperial"){
    stop("CSV writing is only supported for decoupled Imperial sampling. Check human_ode parameter.")
  }

  # check x0
  base_x0(x0)

  # check/setup step function
  stepFun <- base_stepFunc_decoupled(
    sampler = sampler, S = S, hazards = hazards, Sout = Sout,
    dt_stoch = dt_stoch, method = method, inf_labels = inf_labels,
    human_ode=human_ode, ... )

  # check/setup simulation time
  simTimes <- base_time(tt = tmax, dt = dt)

  # check/organize events and x0
  events <- base_events(x0 = x0, events = events, dt = dt)

  # pass everything down to base function - writes CSVs
  sim_trajectory_base_CSV_decoupled(
    x0 = switch(EXPR = sampler, "ode-decoupled" = x0, round(x0)),
    h0 = h0, SPN_P = SPN_P, theta = theta, times = simTimes,
    stepFun = stepFun, events = events, batch = batch, Sout = Sout,
    verbose = verbose, human_ode = human_ode, cube = cube, folders = folders
  )
}


################################################################################
# Base simulation trajectory
################################################################################

#' Simulate Trajectory From one SPN Model using Imperial Malaria Model
#'
#' This is an internal function to \code{\link{sim_trajectory_CSV_decoupled}}. It does the
#' actual sampling once all of the functions have been checked and setup.
#'
#' @param x0 the initial marking of the SPN (initial state)
#' @param h0 initial human state distribution
#' @param SPN_P stochastic Petri Net, places
#' @param theta parameters
#' @param human_ode ode function used for human states
#' @param cube inheritance cube
#' @param times sequence of sampling times
#' @param stepFun a sampling function
#' @param folders vector of folders to write output
#' @param events a \code{data.frame} of events (uses the same format as required
#' in package \code{deSolve} for consistency, see \code{\link[deSolve]{events}}
#' for more information)
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param Sout an optional matrix to track event firings
#' @param verbose print a progress bar?
#'
#' @return no return, prints .csv files into provided folders
#'
#' @importFrom deSolve ode
#' @importFrom stats rbinom
sim_trajectory_base_CSV_decoupled <- function(
    x0, h0, SPN_P, theta, times, stepFun,
    events = NULL, batch = NULL, Sout = NULL, verbose = TRUE,
    human_ode = "Imperial", cube = NULL, folders = folders) {

  # stuff for later
  xNames <- names(x0)
  nTime <- length(times)
  na <- theta$na # number of age compartments
  clin_inc <- theta$clin_inc
  human_state_labels <- generate_Imperial_human_state_labels(na)

  # make index list for all possible output of mosquito epi states
  # humans will be handled differently
  #  stage was hardcoded in Agastya's version, and it was the only change from
  #  the standard function, so removing the duplication and putting
  #  hardcoded params in. - JB
  pStuff <- base_print_stuff(stage = c("F", "M"), state_names = xNames)
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
      pbar <- txtProgressBar(min = 1, max = nTime, style = 3)
      cat(" --- begin simulation --- \n")
    }

    if(human_ode == "Imperial"){
      # initialize human trace list
      names <- as.character(seq(1, 15))
      human_trace <- vector("list", length(names))
      names(human_trace) <- names

      human_trace[["1"]] <- h0
    }

    # reset at the beginning of every simulation
    state <- list("x" = NULL, "o" = NULL, "h" = NULL)
    state$x <- x0
    state$h <- h0
    eventIdx <- 1L
    batchIdx <- 1L


    # Initialize files
    fileCons <- setNames(object = vector(mode = "list", length = lenPS + 1),
                         nm = c(stage, "H"))

    # mosquito files
    for(curS in 1:lenPS){
      # open file connection
      fileCons[[curS]] <- file(description = paste0(folders[num_reps], .Platform$file.sep, stage[curS], ".csv"),
                               open = "wt")
      # write header
      writeLines(text = paste0(c("Time", xNames[pList[[curS]]]), collapse = ","),
                 con = fileCons[[curS]], sep = "\n")
      # write T0
      writeLines(text = paste0(c(0, x0[pList[[curS]]]), collapse = ","),
                 con = fileCons[[curS]], sep = "\n")
    }

    # human file
    # open file connection
    fileCons[["H"]] <- file(description = paste0(folders[num_reps], .Platform$file.sep, "H", ".csv"),
                            open = 'wt')
    # write header
    writeLines(text = paste0(c("Time", human_state_labels), collapse = ","),
               con = fileCons[["H"]], sep = "\n")
    # write T0
    writeLines(text = paste0(c(0, as.vector(h0), clin_inc), collapse = ","),
               con = fileCons[["H"]], sep = "\n")

    # tracking rates/event firing
    if(track){
      event_con <- file(
        description = paste0(folders[num_reps], .Platform$file.sep, "Tracking.csv"),
        open = "wt")
      writeLines(
        text = paste0(c("Time", rownames(Sout)), collapse = ","),
        con = event_con, sep = "\n")
    }

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

      out <- ode(y = state$h, times = c(t0, t1), func = func, parms = theta)

      # update human state matrix
      y <- tail(out, 1)
      num_states <- 10
      # keep clinical incidence separate, it doesn't contribute to ODE dynamics
      # we're just using it for presentation/viz later on
      human_state_matrix <- matrix(y[2:length(y)], ncol = num_states)
      clin_inc_idx <- 10
      clin_inc <- human_state_matrix[, clin_inc_idx]
      human_state_matrix <- human_state_matrix[, -clin_inc_idx]
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

        # fill rest of file with 0s
        #  this makes sure all files are the same size/shape for analysis
        for(nT in i:nTime){
          for(curS in 1:lenPS){
            writeLines(text = paste0(c(times[nT], state$x[pList[[curS]]]), collapse = ","),
                       con = fileCons[[curS]], sep = "\n")
          }
          # human output
          writeLines(
            text = paste0(c(times[nT], rep(0, length(human_state_labels))), collapse = ","),
            con = fileCons[["H"]], sep = ",")
          # tracking rates/event firing
          if(track){
            writeLines(
              text = paste0(c(times[nT], rep(0, nrow(Sout))), collapse = ","),
              con = event_con, sep = "\n")
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
      # mosquitoes
      for(curS in 1:lenPS){
        writeLines(
          text = paste0(c(t1, state$x[pList[[curS]]]), collapse = ","),
          con = fileCons[[curS]], sep = "\n")
      }
      # humans
      writeLines(
        text = paste0(c(t1, as.vector(y[2:length(y)])), collapse = ','),
        con = fileCons[["H"]], sep = "\n")
      # tracking rates/event firing
      if(track){
        writeLines(
          text = paste0(c(t1, state$o), collapse = ","),
          con = event_con, sep = "\n")
      }

      # progress bar
      if(verbose){ setTxtProgressBar(pb = pbar, value = i)}
    } # end sim loop

    # close connections
    # mosquitoes
    for(curS in 1:lenPS){ close(con = fileCons[[curS]]) }
    # humans
    close(con = fileCons[["H"]])
    # tracking
    if(track){ close(con = event_con) }

    if(verbose){
      close(pbar)
      cat(" --- end simulation --- \n")
    }
  } # end repetitions loop

  # no return
}
