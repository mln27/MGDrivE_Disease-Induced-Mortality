################################################################################
#
#   MGDrivE2: SPN structure for a metapopulation network (SEI-SIS epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2019
#
################################################################################

################################################################################
#   Compose transitions to create the metapopulation transition structure
################################################################################

#' Make Transitions (T) For a Network (SEI Mosquitoes - SIS Humans)
#'
#' This function makes the set of transitions (T) for a SPN model of a
#' metapopulation network for simulation of coupled SEI-SIS dynamics. It is the
#' network version of \code{\link{spn_T_epiSIS_node}}.
#'
#' This function takes the places produced from \code{\link{spn_P_epiSIS_network}}
#' and builds all possible transitions between subsets of those places.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SIS}}, \code{\link{equilibrium_SEI_decoupled_mosy}},
#' or \code{\link{equilibrium_Imperial_decoupled}}.
#'
#' While this function produces all structural information related to transitions,
#' hazards are produced by a separate function, \code{\link{spn_hazards}}.
#'
#' For larger networks, this function may take some time to return, please be patient;
#' the Petri Net modeling formalism trades additional computation time at model
#' initialization for faster sampling of trajectories within a simulation.
#'
#' Please note, the movement matrices (\code{h_move} and \code{m_move}) are NOT
#' stochastic matrices, just binary matrices that say if i,j can exchange population.
#' Diagonal elements must be \code{FALSE}, and both matrices are checked for validity; the
#' function will stop with errors if the adjacency matrix specifies illegal movement
#' rules (e.g.; mosquito movement from a "h" node to a "b" node)
#'
#' For examples of using this function, see:
#' \code{vignette("epi-network", package = "MGDrivE2")}
#'
#' @param node_list a character vector specifying what type of nodes to create;
#' (m = a node with only mosquitoes, h = a node with only humans, b = a node with both humans and mosquitoes)
#' @param spn_P set of places produced by \code{\link{spn_P_epiSIS_network}}
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#' @param h_move binary adjacency matrix indicating if movement of humans between nodes is possible or not
#' @param m_move binary adjacency matrix indicating if movement of mosquitoes between nodes is possible or not
#' @param decoupled boolean to denote if using decoupled sims, default is \code{FALSE}
#' @param feqTol tolerance for numerical equality, default is sqrt(.Machine$double.eps)
#'
#' @return a list with two or three elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T), and (if \code{decoupled == TRUE})
#' \code{inf_labels} is the female infection transitions
#'
#' @export
spn_T_epiSIS_network <- function(node_list,spn_P,params,cube,h_move,m_move,decoupled=FALSE,feqTol=1.5e-08){

  # set of places
  u <- spn_P$u

  # genotypes
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # infection stages
  epi_stages <- c("S",paste0("E",as.character(1:params$nEIP)),"I")

  # within node transitions
  T_meta <- vector("list",length(node_list))
  T_index <-  1

  # make all of the transitions in each node
  for(t in 1:length(node_list)){
    if(node_list[t]=="b"){
      T_meta[[t]] <- spn_T_both_epi(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                                    nEIP = params$nEIP,cube = cube,node_id = t,
                                    T_index = T_index, feqTol = feqTol, decoupled = decoupled)
    } else if(node_list[t]=="h"){
      T_meta[[t]] <- spn_T_humans_epi(u = u, node_id = t,T_index = T_index)
    } else if(node_list[t]=="m"){
      T_meta[[t]] <- spn_T_mosy_epi(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                                    nEIP = params$nEIP,cube = cube,node_id = t,
                                    T_index = T_index, feqTol = feqTol)
    } else {
      stop(paste0("error: bad entry in node_list, ",node_list[t]))
    }
  }


  # make mosquito movement events
  # This checks for valid movement, then sets lists for male/female mosquitoes,
  #  and human movement
  if(any(Matrix::diag(h_move))){
    stop("adjacency matrix 'h_move' must have all FALSE elements along its diagonal")
  }
  if(any(Matrix::diag(m_move))){
    stop("adjacency matrix 'm_move' must have all FALSE elements along its diagonal")
  }
  moveList <- spn_T_move_epi(node_list = node_list,u = u,m_move = m_move,
                             h_move = h_move,nG = nG,g = g,nEIP = params$nEIP,
                             epi_stages = epi_stages,
                             h_state = c("S","I"),T_index = T_index)


  # set of transitions with movement appended to it
  v <- c(unlist(x = lapply(X = T_meta, FUN = '[[', 'v'), use.names = FALSE),
         unlist(x = lapply(X = moveList, FUN = lapply, '[[', 'label'),
                use.names = FALSE))

  # all the transitions (and flatten them)
  T_all <- c(list("T_win" = do.call(c,lapply(X = T_meta, FUN = '[[', 'T'))),
            moveList)

  T_all <- unlist(x = T_all, recursive = FALSE, use.names = FALSE)


  # check the set for errors
  labels <- unlist(x = lapply(X = T_all, FUN = '[[', 'label'),use.names = FALSE)
  indices <- unlist(x = lapply(X = T_all, FUN = '[[', 'vix'),use.names = FALSE)
  if(any(labels != v[indices])){
    stop(paste0("error in set of transitions T and transition vector at transition(s): ",
                paste0(v[(labels != v[indices])], collapse = ", ")))
  }


  # return the complete set of transitions {T}
  return(list("T" = T_all,
              "v" = v) )
}


################################################################################
# WITHIN NODE INTERNAL TRANSITION FUNCTIONS
#
# NOTE:
# these functions make the within node transitions for the T set.
# they are not exported from the package for users
# as the package provides support only to generate Petri Nets by
# returning complete sets (P,T); if the below functions were exported
# it would be possible for a user to return an incomplete set of places (P)
# or transitions (T)
#
# the reason we don't use the {T} generator for the one-node epi
# system here is because there would be too much obfuscating logic in the
# individual transition functions to check for null nodes, etc.
#
################################################################################

# move one human
make_transition_move_human_epiSIS <- function(T_index,u,state,origin,dest){

  # safety check
  stopifnot(state %in% c("S","E","I","R"))

  # tokens required/produced
  htoken1_2 <- file.path("H", state, c(origin,dest), fsep = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(htoken1_2, collapse = "->") # name of this t (corresponds to v)

  # requires one human token
  t$s <- match(x = htoken1_2[1], table = u)
  t$s_w <- 1

  # produces one human token
  t$o <- match(x = htoken1_2[2], table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "move_human"

  # return the transition
  return(t)
}

################################################################################
# make movement transitions (T) for everyone SIS
################################################################################

spn_T_move_epi <- function(node_list,u,m_move,h_move,nG,g,nEIP,epi_stages,
                           h_state,T_index){

  # before we set up movement events, we run the checker to make sure we don't make illegal movement:
  # eg; humans can't move to mosquito only locations and vice versa
  check_move_legal(node_list=node_list,move=h_move,checker=h_move_check,hm="human")
  check_move_legal(node_list=node_list,move=m_move,checker=m_move_check,hm="mosquito")

  # make mosquito movement events

  # edges: need to make male and female movement for these edges
  m_edges <- arrayInd(Matrix::which(m_move),.dim=dim(m_move))
  m_edges <- m_edges[order(m_edges[,1]),]

  # make female movement
  female_move <- vector("list",nrow(m_edges)*(nG^2)*(nEIP+2))
  vv <- 1

  # iterate over edges
  for(e in 1:nrow(m_edges)){
    # ... and over infection status
    for(s in 1:(nEIP+2)){
      # ... and over female genotypes
      for(f in 1:nG){
        # ... and over male mate genotypes
        for(m in 1:nG){
          female_move[[vv]] <- make_transition_move_female(T_index=T_index,
                                                               u=u,f_gen=g[f],
                                                               m_gen=g[m],
                                                               inf=epi_stages[s],
                                                               origin=m_edges[e,1],
                                                               dest=m_edges[e,2])
          T_index <- T_index + 1
          vv <- vv + 1
        }
      }
    }
  }

  # make male movement
  male_move <- vector("list",nrow(m_edges)*nG)
  vv <- 1

  # iterate over edges
  for(e in 1:nrow(m_edges)){
    # ... and over male genotypes
    for(m in 1:nG){
      male_move[[vv]] <- make_transition_move_male(T_index = T_index,u = u,
                                                   m_gen = g[m],
                                                   origin = m_edges[e,1],
                                                   dest = m_edges[e,2])
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # make human movement events

  # edges: need to make male and female movement for these edges
  h_edges <- arrayInd(Matrix::which(h_move),.dim=dim(h_move))
  h_edges <- h_edges[order(h_edges[,1]),]

  # make human movement
  human_move <- vector("list",nrow(m_edges)*length(h_state))
  vv <- 1

  # iterate over edges
  for(e in 1:nrow(h_edges)){
    for(state in h_state){
      human_move[[vv]] <- make_transition_move_human_epiSIS(T_index=T_index,u=u,
                                                          state=state,
                                                          origin=h_edges[e,1],
                                                          dest=h_edges[e,2])
      T_index <- T_index + 1
      vv <- vv + 1
    } # end loop over state
  } # end movement loop


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return named list
  return(list("female_move" = female_move,
              "male_move" = male_move,
              "human_move" = human_move) )
}

################################################################################
# make base internal transitions (T) for other functions below
################################################################################

# base mosy functions
#  this is used in spn_T_mosy_epi and spn_T_both_epi
# This function takes a bunch of shit
#  Returns a list of lists with all of the base actions
base_T_mosy_epi <- function(u,nE,nL,nP,nEIP,cube,node_id,T_index,epi_stages,feqTol){

  # genetic states
  g <- cube$genotypesID
  nG <- cube$genotypesN

  # base mosquito stuff
  trans <- base_T_mosy(u=u,nE=nE,nL=nL,nP=nP,nG=nG,g=g,node_id=node_id,T_index=T_index)

  # empty list to put the transitions in (X_tt is the set of transitions in this subset of the total T)
  epi_len <- length(epi_stages)
  ovi_tt_len <- sum((cube$tau * cube$ih) > 0) * epi_len
  ovi_tt <- vector("list", ovi_tt_len)
  vv <- 1

  # OVIPOSITION

  # make the transitions
  ovi_dims <- dim(cube$ih)
  for(i in 1:ovi_dims[1]){
    for(j in 1:ovi_dims[2]){
      for(k in 1:ovi_dims[3]){
        # only make valid events (based on tau)
        if(!fequal(x = cube$ih[i,j,k], y = 0, tol = feqTol) &&
           !fequal(x = cube$tau[i,j,k], y = 0, tol = feqTol)){
          # loop over epi stages
          for(l in 1:epi_len){
            ovi_tt[[vv]] <- make_transition_ovi_epi(T_index,u=u,f_gen=g[i],
                                                    m_gen=g[j],o_gen=g[k],
                                                    inf=epi_stages[l],node=node_id)
            T_index <- T_index + 1
            vv <- vv + 1
          }
        }
      }
    }
  }

  # check that both checks return the same thing
  if(ovi_tt_len != (vv-1)){
    stop(paste0("Some values of tau or ih are less than ", feqTol,
           ".\n\tIf this is desired, please reduce PARAM"))
  }


  # PUPAE TRANSITIONS

  # make pupae -> female emergence
  pupae_2female_tt <- vector("list",nG^2)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      pupae_2female_tt[[vv]] <- make_transition_pupae_emerge_f_epi(T_index,
                                                                   u=u,p_gen=g[j_f],
                                                                   m_gen=g[j_m],
                                                                   nP=nP,node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # FEMALE TRANSITIONS

  # make female mortality
  female_mort_tt <- vector("list",(nG^2)*length(epi_stages))
  vv <- 1

  for(j_eip in 1:length(epi_stages)){
    for(j_f in 1:nG){
      for(j_m in 1:nG){
        female_mort_tt[[vv]] <- make_transition_female_mort_epi(T_index,u=u,f_gen=g[j_f],
                                                                m_gen=g[j_m],
                                                                inf=epi_stages[j_eip],
                                                                node=node_id)
        T_index <- T_index + 1
        vv <- vv + 1
      }
    }
  }

  # UNMATED FEMALE MATING
  unmated_mate_tt <- vector("list",nG)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      unmated_mate_tt[[vv]] = make_transition_female_unmated_mate(T_index,u=u,
                                                                  f_gen=g[j_f],
                                                                  m_gen=g[j_m],
                                                                  node=node_id,
                                                                  epi=TRUE)
      T_index = T_index + 1
      vv = vv + 1
    }
  }

  # no female infection (no humans here!)

  # make female eip (incubation)
  female_eip_tt <- vector("list",0)

  if(nEIP > 1){

    female_eip_tt <- vector("list",(nG^2)*(nEIP-1))
    vv <- 1

    for(j_eip in 1:(nEIP-1)){
      for(j_f in 1:nG){
        for(j_m in 1:nG){
          female_eip_tt[[vv]] <- make_transition_female_eip_epi(T_index,u=u,
                                                                     f_gen=g[j_f],
                                                                     m_gen=g[j_m],
                                                                     inc1=j_eip,
                                                                     inc2=j_eip+1,
                                                                     node=node_id)
          T_index <- T_index + 1
          vv <- vv + 1
        }
      }
    }
  }

  # make female infectiousness (clear eip)
  female_inc_tt <- vector("list",nG^2)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      female_inc_tt[[vv]] <- make_transition_female_inc_epi(T_index,u=u,
                                                                 f_gen=g[j_f],
                                                                 m_gen=g[j_m],
                                                                 nEIP=nEIP,node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }


  # push T_index back up
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return as a list for further processing
  return(c(trans,
           list("oviposit" = ovi_tt,
                "pupae_2female" = pupae_2female_tt,
                "female_mort" = female_mort_tt,
                "unmated_mate" = unmated_mate_tt,
                "female_eip" = female_eip_tt,
                "female_inc" = female_inc_tt)
           ) )

} # end base mosy func

# base human functions
#  this is used in spn_T_humans_epi and spn_T_both_epi
# This function takes a bunch of shit
#  Returns a named list with all the base transitions
base_T_humans_epi <- function(u,node_id,T_index){

  # HUMAN TRANSITIONS
  human_tt <- vector("list",5)
  vv <- 1

  human_tt[[vv]] <- make_transition_human_birth_S_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_birth_I_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_death_S_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_death_I_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_rec_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  return(human_tt)
}

################################################################################
# make internal transitions (T) for mosquito-only node
################################################################################

# u: set of places
# if this is for node 1, T_index = 1, otherwise its max(node[i-1].T_index)+1
# node_id: what node is this
spn_T_mosy_epi <- function(u,nE,nL,nP,nEIP,cube,node_id,T_index,feqTol){

  epi_stages <- c("S",paste0("E",as.character(1:nEIP)),"I")

  # make oviposition transitions (events)
  t <- base_T_mosy_epi(u =u,nE= nE,nL= nL,nP = nP,nEIP= nEIP,cube = cube,
                       node_id=node_id,T_index = T_index,epi_stages = epi_stages,
                       feqTol = feqTol)

  # transitions (v)
  v <- unlist(x = lapply(X = t, FUN = lapply, '[[', 'label'), use.names = FALSE)

  # one long vector
  t <- unlist(x = t, recursive = FALSE, use.names = FALSE)


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return the set of transitions and the vector (v)
  return(list("T" = t,
              "v" = v) )
}

################################################################################
# make internal transitions (T) for human-only node
################################################################################

# u: set of places for the mosquito-only node
# if this is for node 1, T_index = 1, otherwise its max(node[i-1].T_index)+1
# node_id: id of the node these transitions are for
spn_T_humans_epi <- function(u,node_id,T_index){

  # HUMAN TRANSITIONS
  human_tt <- base_T_humans_epi(u = u,node_id = node_id ,T_index = T_index)

  # transitions (v)
  v <- unlist(x = lapply(X = human_tt, FUN = '[[', 'label'), use.names = FALSE)


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return the set of transitions and the vector (v)
  return(list("T" = human_tt,
              "v" = v) )
}


################################################################################
# make internal transitions (T) for human & mosquito node
################################################################################

# u: set of places for the mosquito-only node
# if this is for node 1, T_index = 1, otherwise its max(node[i-1].T_index)+1
spn_T_both_epi <- function(u,nE,nL,nP,nEIP,cube,node_id,T_index,feqTol,decoupled){

  # genetic states
  g <- cube$genotypesID
  nG <- cube$genotypesN

  epi_stages <- c("S",paste0("E",as.character(1:nEIP)),"I")

  # MOSQUITO TRANSITIONS
  # base things in a list
  base_mos <- base_T_mosy_epi(u = u, nE = nE, nL = nL, nP = nP, nEIP = nEIP,
                               cube = cube, node_id = node_id, T_index = T_index,
                               epi_stages = epi_stages, feqTol = feqTol)

  # make female infection
  female_inf_tt <- vector("list",nG^2)
  vv <- 1
  inf_func <- if(!decoupled) make_transition_female_inf_epi else make_transition_female_inf_epi_decoupled

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      female_inf_tt[[vv]] <- inf_func(T_index,u=u,f_gen=g[j_f],m_gen=g[j_m],node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # check for normal or decoupled sims
  if(!decoupled){
    # running normal sims, full Petri Net

    # HUMAN TRANSITIONS
    human_tt <- vector("list",nG^2)
    vv <- 1

    # add epi transitions not in base
    for(j_f in 1:nG){
      for(j_m in 1:nG){
        human_tt[[vv]] <- make_transition_human_inf_epiSIS(T_index,u=u,f_gen=g[j_f],
                                                           m_gen=g[j_m],node=node_id)
        T_index <- T_index + 1
        vv <- vv + 1
      }
    }

    # combine with base transitions
    human_tt <- c(human_tt,
                  base_T_humans_epi(u = u, node_id = node_id, T_index = T_index))

    # the set of transitions
    t <- c(base_mos,
           "female_inf" = list(female_inf_tt),
           "human" = list(human_tt))

  } else {
    # decoupled sims, for Imperial model

    # the set of transitions
    t <- c(base_mos,
           "female_inf" = list(female_inf_tt))

    # keep track of the female infection transitions
    # since these rely on human states, we will need them later
    # when evaluating hazards
    inf_labels <- unlist(sapply(female_inf_tt, function(inf) { return(inf$label) }))
  }


  # transitions (v)
  v <- unlist(x = lapply(X = t, FUN = lapply, '[[', 'label'), use.names = FALSE)

  # one long vector
  t <- unlist(x = t, recursive = FALSE, use.names = FALSE)

  # setup return
  retList <- if(!decoupled){
    list("T" = t,"v" = v)
  } else {
    list("T" = t,"v" = v,"inf_labels" = inf_labels)
  }


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return the set of transitions and the vector (v)
  return(retList)
}
