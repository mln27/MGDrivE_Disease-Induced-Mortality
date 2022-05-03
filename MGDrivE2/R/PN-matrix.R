################################################################################
#
#   MGDrivE2: generate sparse matrices from PN
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#   Jared Bennett (jared_bennett@berkeley.edu)
#   May 2022
#
################################################################################

################################################################################
# make the Pre matrix (v by u)
################################################################################

#' Make Pre Matrix For a Petri Net
#'
#' Generate the Pre (|v| by |u|) matrix for the SPN. This gives the edges from P
#' to T (input arcs) in the bipartite network.
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
#'
#' @param spn_P set of places (P) (see details)
#' @param spn_T set of transitions (T) (see details)
#'
#' @return a matrix of type \code{\link[Matrix]{dgCMatrix-class}}
#'
#' @importFrom Matrix sparseMatrix
#'
spn_Pre <- function(spn_P,spn_T){
  # These updates owe a lot to Sean and his vectorization in the dev branch
  # Some of those ideas came from me, the implementation and finishing are him
  # The updates below modify the core logic but maintain much of his work

  # these will be used at the end
  u <- spn_P$u # dimension of the places
  v <- spn_T$v # dimension of the transitions

  # indexes, input arcs, and input weights
  ix <- lapply(X = spn_T$T, FUN = "[[", "vix")
  s <- lapply(X = spn_T$T, FUN = "[[", "s")
  s_w <- lapply(X = spn_T$T, FUN = "[[", "s_w")

  # ignore transitions that are always on (if any)
  #  This replaces the heart of the loop from the original function
  #  It uses lapply because it was faster (empirical testing)
  #  The logic is inverted, because it is faster (more/earlier escape places)
  #   and we don't have to invert it again afterward - we're not really interested
  #   in the "FALSE" indexes, we want to keep the "TRUE" ones.
  #  unlist to get vector
  #  which so we can do index-based subsetting - way faster
  not_null_idx <- which(unlist(lapply(X = seq_len(length(s)),
                                      FUN = function(i){any(!is.nan(s[[i]]), !is.nan(s_w[[i]]))}
                                      )))

  # keep things
  #  notice index-based subsetting
  ix <- ix[not_null_idx]
  s <- s[not_null_idx]
  s_w <- s_w[not_null_idx]

  # for transitions with multiple input arcs
  s_len <- lengths(x = s, use.names = FALSE)

  # replicate the elements of the ix (transition index) that number of times
  #  we replicate everything because it is faster than handling pieces in R
  #  This also turns "ix" into an integer vector, necessary for the matrix
  ix <- rep.int(x = as.integer(unlist(ix)), times = s_len)

  # build sparse pre-matrix
  #  this is transposed from the original - from v X u matrix to a u X v, because
  #  that's the final shape we want
  Pre <- sparseMatrix(
    i = as.integer(unlist(s)),
    j = ix,
    x = as.integer(unlist(s_w)),
    dims = c(length(u), length(v)),
    dimnames = list(u, v)
  )

  return(Pre)
}


################################################################################
# make the Post matrix (v by u)
################################################################################

#' Make Post Matrix For a Petri Net
#'
#' Generate the Post (|v| by |u|) matrix for the SPN. This gives the edges from
#' T to P (output arcs) in the bipartite network.
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
#'
#' @param spn_P set of places (P) (see details)
#' @param spn_T set of transitions (T) (see details)
#'
#' @return a matrix of type \code{\link[Matrix]{dgCMatrix-class}}
#'
spn_Post <- function(spn_P,spn_T){
  # These updates owe a lot to Sean and his vectorization in the dev branch
  # Some of those ideas came from me, the implementation and finishing are him
  # The updates below modify the core logic but maintain much of his work

  # these will be used at the end
  u <- spn_P$u # dimension of the places
  v <- spn_T$v # dimension of the transitions

  # indexes, output arcs, and output weights
  ix <- lapply(X = spn_T$T, FUN = "[[", "vix")
  o <- lapply(X = spn_T$T, FUN = "[[", "o")
  o_w <- lapply(X = spn_T$T, FUN = "[[", "o_w")

  # ignore transitions that are always on (if any)
  #  This replaces the heart of the loop from the original function
  #  It uses lapply because it was faster (empirical testing)
  #  The logic is inverted, because it is faster (more/earlier escape places)
  #   and we don't have to invert it again afterward - we're not really interested
  #   in the "FALSE" indexes, we want to keep the "TRUE" ones.
  #  unlist to get vector
  #  which so we can do index-based subsetting - way faster
  not_null_idx <- which(unlist(lapply(X = seq_len(length(o)),
                                      FUN = function(i){any(!is.nan(o[[i]]), !is.nan(o_w[[i]]))}
  )))

  # keep things
  #  notice index-based subsetting
  ix <- ix[not_null_idx]
  o <- o[not_null_idx]
  o_w <- o_w[not_null_idx]

  # for transitions with multiple output arcs
  o_len <- lengths(x = o, use.names = FALSE)

  # replicate the elements of the ix (transition index) that number of times
  #  we replicate everything because it is faster than handling pieces in R
  #  This also turns "ix" into an integer vector, necessary for the matrix
  ix <- rep.int(x = as.integer(unlist(ix)), times = o_len)

  # build sparse post-matrix
  #  this is transposed from the original - from v X u matrix to a u X v, because
  #  that's the final shape we want
  Post <- sparseMatrix(
    i = as.integer(unlist(o)),
    j = ix,
    x = as.integer(unlist(o_w)),
    dims = c(length(u), length(v)),
    dimnames = list(u, v)
  )

  return(Post)
}


################################################################################
# make the other matrices that are useful
# A: reaction matrix; v by u
# S: stoichiometry matrix; u by v
################################################################################

#' Make stoichiometry Matrix For a Petri Net
#'
#' Generate the stoichiometry (|u| by |v|) matrix for the SPN.
#' Each column gives the net effect of that transition firing upon the state
#' space of the model. Internally, this creates a Pre (\code{\link{spn_Pre}}) and
#' Post (\code{\link{spn_Post}}) matrix, and then calculates the final stoichiometry.
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
#'
#' @param spn_P set of places (P) (see details)
#' @param spn_T set of transitions (T) (see details)
#'
#' @importFrom Matrix drop0
#'
#' @export
spn_S <- function(spn_P,spn_T){

  # create pre matrix
  Pre <- spn_Pre(spn_P = spn_P,spn_T = spn_T)

  # create post matrix
  Post <- spn_Post(spn_P = spn_P,spn_T = spn_T)

  # calculate difference stoichiometry
  #  drop zeros
  #  return (u X v) matrix
  # A matrix
  return(drop0(x = Post - Pre))
}
