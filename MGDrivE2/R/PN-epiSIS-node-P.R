################################################################################
#
#   MGDrivE2: SPN structure for a single node (SEI-SIS epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################


################################################################################
# make the places (P) of the SPN
################################################################################

#' Make Places (P) For a Node (SEI Mosquitoes - SIS Humans)
#'
#' This function makes the set of places (P) for a SPN. It is used alone if our model
#' is a single-node metapopulation for mosquito SEI and human SIS dynamics;
#' otherwise it is used as part of other functions to make SPN models with
#' larger state spaces (metapopulation models, see \code{\link{spn_P_epiSIS_network}}).
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SIS}}, \code{\link{equilibrium_SEI_decoupled_mosy}}
#' or \code{\link{equilibrium_Imperial_decoupled}}.
#'
#' If it is used with the decoupled simulations, we only
#' use the SPN framework for the mosquito component of the model. The human component
#' will be handled separately in the sampler, and will be formulated as an ODE.
#' This is used by both SIS-decoupled and Imperial transmission models
#'
#' For examples of using this function, see:
#' \code{vignette("epi-node", package = "MGDrivE2")},
#' \code{vignette("epi-node-decoupled", package = "MGDrivE2")} or
#' \code{vignette("epi-node-imperial", package = "MGDrivE2")}.
#'
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#' @param decoupled boolean to denote if using decoupled sims, default is \code{FALSE}
#'
#' @return a list with two elements: \code{ix} contains labeled indices of the places
#' by life stage, \code{u} is the character vector of places (P)
#'
#' @export
spn_P_epiSIS_node <- function(params,cube,decoupled=FALSE){

  # checks
  nE <- params$nE
  nL <- params$nL
  nP <- params$nP
  nEIP <- params$nEIP
  if(nE < 2 || nL < 2 || nP < 2 || nEIP < 2){
    warning(paste0("A shape parameter ('nE', 'nL', 'nP', 'nEIP') of 1 implies ",
                   "exponentially distributed dwell times in that compartment."))
  }

  # genetic information
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # setup places
  eggs <- file.path("E",1:nE,"_",rep(g, each = nE), fsep = "")

  larvae <- file.path("L",1:nL,"_",rep(g, each = nL), fsep = "")

  pupae <- file.path("P",1:nP,"_",rep(g, each = nP), fsep = "")

  females_unmated <- file.path("U",g, fsep = "_")

  # create female/male geno combinations, then replicate that for each EIP
  # tack on infection status, replicated for the number of genotypes.
  stages <- c("S",paste0("E",as.character(1:nEIP)),"I")
  females <- file.path("F", rep(x = rep(x = g, each = nG), times = nEIP+2),
                      rep(x = g, times = nG*(nEIP+2)),
                      rep(x = stages, each = nG^2),
                      fsep = "_")

  males <- file.path("M",g, fsep = "_")

  # indices of states
  ix <- list()
  ix$egg <- matrix(seq_along(eggs),nrow = nE,byrow = FALSE,dimnames = list(1:nE,g))
  ix$larvae <- matrix(data = seq_along(larvae) + nG*nE,nrow = nL,byrow = FALSE,dimnames = list(1:nL,g))
  ix$pupae <- matrix(data = seq_along(pupae) + nG*(nE + nL),nrow = nP,byrow = FALSE,dimnames = list(1:nP,g))
  ix$females_unmated <- setNames(object = seq_along(females_unmated) + nG*(nE + nL + nP), nm = g)

  # For the females, we have to permute because the array function doesn't have
  # something like "byrow" for matrices. This makes sure the "rows" are the female
  # genotype, cols are male mate, and z-axis is stage of incubation.
  # "resize=TRUE" is a cheat to keep dim names, doesn't cause issues in this instance
  ix$females <- aperm(a = array(data = seq_along(females) + nG*(nE + nL + nP + 1),dim = c(nG,nG,nEIP+2),
                  dimnames = list(g,g,stages)),
        perm = c(2,1,3),resize = TRUE)

  ix$males <- setNames(object = seq_along(males) + nG*(nE+nL+nP+nG*(nEIP+2) + 1), nm = g)

  # places (u)
  u <- c(eggs,larvae,pupae,females_unmated,females,males)

  # check for decoupled setup
  if(!decoupled){
    # we're doing normal sims, calculate and add humans
    # human places
    humans <- c("H_S","H_I")
    # human indices
    ix$humans <- setNames(object = seq_along(humans) + nG*(nE+nL+nP+nG*(nEIP+2) + 2),nm = humans)
    # add humans to places
    u <- c(u,humans)
  }

  # return list of places
  #  make ix a list to match network version
  return(list("ix" = list(ix),
              "u" = u) )
}
