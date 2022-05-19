###############################################################################
#
#   MGDrivE2: hazard functions for a metapopulation network (SEI-SEIR epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
###############################################################################

###############################################################################
#   Make the SPN hazards
###############################################################################

#' Make Hazards (Lambda) For a MGDrivE2: Node and Network Simulations
#'
#' Using the structural (topological) SPN model as well as parameters in the
#' \code{cube} and \code{params} objects, generate a list (of length |v|) of
#' hazards, each implemented as a function closure.
#'
#' If these hazards will be used in a continuous approximation algorithm, such as
#' an ODE method (\code{\link{step_ODE}}) or Gillespie's Direct Method
#' (\code{\link{step_DM}}), it is recommended to use \code{exact=FALSE}. If the
#' hazards will be used in an integer state space method, such as tau-leaping
#' (\code{\link{step_PTS}}) or Chemical Langevin (\code{\link{step_CLE}}) methods,
#' it is recommended to use \code{exact=TRUE}.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The set of transitions (\code{spn_T}) is generated from one of the following:
#' \code{\link{spn_T_lifecycle_node}}, \code{\link{spn_T_lifecycle_network}},
#' \code{\link{spn_T_epiSIS_node}}, \code{\link{spn_T_epiSIS_network}},
#' \code{\link{spn_T_epiSEIR_node}}, \code{\link{spn_T_epiSEIR_network}}.
#'
#' The \code{params} objected is generated from either \code{\link{equilibrium_lifeycle}}
#' or \code{\link{equilibrium_SEI_SIS}}; it is the "params" object in the return
#' list. The equilibrium function used must match the \code{type} parameter.
#'
#' The \code{type} parameter indicates what type of simulation is being run. It
#' is one of: "life", "SIS", or "SEIR". This must match the \code{params} object
#' supplied.
#'
#' Use of this function is demonstrated in many vignettes, \code{browseVignettes(package = "MGDrivE2")}
#'
#' @param spn_P the set of places (P) (see details)
#' @param spn_T the set of transitions (T) (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#' @param params a named list of parameters (see details)
#' @param type string indicating type of hazards, one of; "life", "SIS", or "SEIR"
#' @param log_dd if \code{TRUE}, use logistic (carrying capacity) density dependent hazards, if \code{FALSE} use Lotka-Volterra density dependent hazards for larval mortality
#' @param exact boolean, make exact (integer input) hazards? Default is TRUE
#' @param tol if \code{exact=FALSE}, the value of hazard below which it is clipped to 0
#' @param verbose display a progress bar when making hazards?
#'
#' @return list of length 3: \code{hazards} is a list of named closures for every
#' state transition in the model, \code{flag} is a boolean indicating exact or approximate,
#' \code{approxS} is a vector of state dependencies used in the first-order
#' sparsity-aware samplers
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
spn_hazards <- function(spn_P,spn_T,cube,params,type="life",
                        log_dd=TRUE,exact=TRUE,tol=1e-12,verbose=TRUE){

  # approx vs exact
  check_approx(tol = tol, exact = exact)

  # check parameters are properly specified
  if(type == "life"){
    check_params_life(params = params)
  } else if(type == "SIS"){
    check_params_sis(params = params)
  } else if(type == "SEIR"){
    check_params_SEIR(params = params)
  } else {
    stop("type of hazards improperly specified, please provide: life, SIS, SEIR")
  }

  # check type/validity of density dependence
  check_dd(log_dd = log_dd, params = params)

  # check omega
  #  is.vector() is true for vectors and lists
  #  is.list() is not true for vectors
  #  https://stackoverflow.com/questions/19501186/how-to-test-if-object-is-a-vector
  if(!is.list(cube$omega)){
    # basic omega, tweak for infection status
    cube$omega <- list("S"=cube$omega,
                       "E"=cube$omega,
                       "I"=cube$omega)
  } else {
    # omega already set for infection status
    #  check for lengths and names
    #  these are equivalent to the checks performed on omega in the cube

    # lengths
    if(!all(lengths(cube$omega)==cube$genotypesN)){
      stop("length(s) of omega do not match the number of genotypes in the cube;
           please check the length(s) in omega")
    }

    # names
    if(!all(vapply(X = 1:length(cube$omega),
                   FUN = function(x){
                     all(names(cube$omega[[x]]) %in% cube$genotypesID)},
                   FUN.VALUE = logical(length = 1L))
    )){
      stop("genotype(s) in omega do not match genotypes in cube;
           please check names of omega")
    }

    # values must be >0
    if(any(unlist(x = cube$omega, use.names = FALSE) < 0)){
      stop("omega values must be positive [X>0]")
    }
  } # end omega reshape and check


  # transitions and places
  v <- spn_T$v
  u <- spn_P$u

  n <- length(v)
  if(verbose){
    pb <- txtProgressBar(min = 1,max = n,style = 3)
    pp <- 1

    cat(" --- generating hazard functions for SPN --- \n")
  }

  # the hazard functions
  h <- setNames(object = vector("list",n),nm = v)


  # setup list of functions

  # basic life cycle
  classNames <- c("oviposit","egg_adv","egg_mort","larvae_adv","pupae_adv",
                  "pupae_mort","pupae_2m","male_mort","female_mort",
                  "female_unmated_mort" )
  funcs <- setNames(object = c(make_oviposit_haz,make_egg_adv_haz,make_egg_mort_haz,
                               make_larvae_adv_haz,make_pupae_adv_haz,make_pupae_mort_haz,
                               make_pupae_2male_haz,make_male_mort_haz,make_female_mort_haz,
                               make_female_mort_haz),
                    nm = classNames)

  # larvae
  lFuncs <- c(make_larvae_mort_haz_log,make_larvae_mort_haz_lk)

  # mating: these need the indices of larvae in the node
  matingNames <- c("pupae_2f","pupae_2unmated","female_unmated_mate")
  matingFuncs <- c("pupae_2f" = make_pupae_2female_haz,
                   "pupae_2unmated" = make_pupae_2unmated_haz,
                   "female_unmated_mate" = make_unmated_2female_haz)

  # infection (M->H or H->M): these need the indices of humans in the node
  infNames <- c("female_inf","H_infection")
  infFuncs <- setNames(object = c(make_female_inf_epi_haz,make_human_inf_sis_haz),
                      nm = infNames)

  # epidemiological hazards
  eipHNames <- c("female_eip","female_inc","H_birth","H_mort","H_recovery",
                  "move_male","move_female","move_human","H_latent")
  eipHFuncs <- setNames(object = c(make_female_eip_epi_haz,make_female_eip_epi_haz,
                                  make_human_birth_sis_haz,make_human_death_sis_haz,
                                  make_human_rec_sis_haz,make_move_male_haz,
                                  make_move_female_haz,make_human_move_sis_haz,
                                  make_human_latent_seir_haz),
                        nm = eipHNames)


  # make the hazards
  for(cT in 1:n){

    type <- spn_T$T[[cT]]$class

    # make the correct type of hazard
    if(type %in% classNames){
      h[[cT]] <- funcs[[type]](trans = spn_T$T[[cT]],u = u,cube = cube,params = params,exact = exact,tol = tol)
    } else if(type == "larvae_mort"){
      # get larvae indices
      #  don't we know the pattern of this?
      #  ie, L#_geno_node, so we always want the 3rd split
      node <- suppressWarnings(as.integer(tail(x = strsplit(x = u[spn_T$T[[cT]]$s],split = "_",fixed = TRUE)[[1]],
                              n = 1)))

      # safety for 1-node "network"
      if(is.na(node)){node <- 1}

      l_ix <- as.vector(spn_P$ix[[node]]$larvae)

      # need true to be 1, false to be 2, so invert bool and cast as double
      # FAlse == 0, so !TRUE + 1 = 1
      h[[cT]] <- lFuncs[[(!log_dd) + 1]](trans = spn_T$T[[cT]],u = u,l_ix = l_ix,node=node,
                                        cube = cube,params = params,exact = exact,tol = tol)
    } else if(type %in% matingNames){
      # get male indices
      #  don't we know the pattern of this?
      #  ie, P#_geno_node, so we always want the 3rd split
      node <- suppressWarnings(as.integer(tail(x = strsplit(x = u[spn_T$T[[cT]]$s[1]],split = "_",fixed = TRUE)[[1]],
                              n = 1)))

      # safety for 1-node "network"
      if(is.na(node)){node <- 1}

      m_ix <- spn_P$ix[[node]]$males

      # make the hazard
      h[[cT]] <- matingFuncs[[type]](trans = spn_T$T[[cT]],u = u,m_ix = m_ix,
                                       cube = cube,params = params,exact = exact,
                                       tol = tol)
    } else if(type %in% infNames){
      # get human indices
      #  don't we know the pattern of this?
      node <- suppressWarnings(as.integer(tail(x = strsplit(x = u[spn_T$T[[cT]]$s[1]],split = "_",fixed = TRUE)[[1]],
                              n = 1)))

      # safety for 1-node "network"
      if(is.na(node)){node <- 1}

      h_ix <- spn_P$ix[[node]]$humans

      # make the hazard
      h[[cT]] <- infFuncs[[type]](trans = spn_T$T[[cT]],u = u,h_ix = h_ix,cube = cube,
                                 params = params,exact = exact,tol = tol)
    } else if(type %in% eipHNames){
      h[[cT]] <- eipHFuncs[[type]](trans = spn_T$T[[cT]],u = u,params = params,exact = exact,tol = tol)
    }  else {
      stop(paste0("error in making hazard function for unknown class type: ",type))
    }


    if(verbose){setTxtProgressBar(pb,cT)}
  } # end loop

  if(verbose){
    close(pb)
    cat(" --- done generating hazard functions for SPN --- \n")
  }

  # Sparsity Indices
  #  Each transition depends on the value of a specific state in the simulation
  #  None of our hazards create something from nothing - so empty states mean 0 hazards
  #  If we keep track of those state dependencies, and then check them for value,
  #  we could skip any state that leads to a 0 hazard without needing to calculate
  #  the hazard. This basically skips any hazards that depend on states with no counts
  #  at a given time-step.
  #  This will fail if we have hazards that create value from nothing
  #   (currently, none of these exist in the sim)
  #  Some hazards depend on multiple input states - namely, the mating functions,
  #   which rely on female and male states:
  #    "pupae_2f" = make_pupae_2female_haz
  #    "female_unmated_mate" = make_unmated_2female_haz
  #   However, if either state is 0, then
  #   the hazard is 0, so by checking only 1 of the 2 states, we are either right
  #   (value of state is 0, and we skip evaluation of the hazard), or the checked
  #   state has value but the unchecked state does not, in which case we needlessly
  #   evaluate a hazard that returns zero. Thus, no statistical issues, but some
  #   loss of performance.
  #  This pulls out the first value of every hazard's state dependence
  approxS <- vapply(X = 1:n, FUN = function(x){spn_T$T[[x]]$s[1]}, FUN.VALUE = integer(length = 1))


  return(list("hazards"=h,"flag"=exact,"approxS"=approxS))
}
