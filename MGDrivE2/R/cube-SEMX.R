###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE2: Mosquito Gene Drive Explorer V2 with Disease-Induced Mortality
#   ReMEDE (REpeat Mediated Excision of a Drive Element) X-Linked
#   Jared Bennett
#   jared_bennett@berkeley.edu
#   March 2022
#
###############################################################################

#' Inheritance Cube: X-Linked ReMEDE (REpeat Mediated Excision of a Drive Element)
#'
#' The ReMEDE system, put forth by \href{GETLINKWHENEXISTS!!!}{Chennuri and Myles},
#' is a realized application of \href{https://doi.org/10.1098/rstb.2019.0804}{biodegradable gene drives}
#' with application of \href{https://doi.org/10.1016/j.celrep.2020.107841}{Small-Molecule Control}
#' for the removal of drive elements through induceable SSA.
#'
#' This construct has a
#' single, X-Linked target site that carries gRNA, Cas9, an induceable endonuclease,
#' and is flanked by direct repeats. There are 8 possible alleles:
#' \itemize{
#'  \item X: Wild-type X chromosome
#'  \item G: Active GD (gene drive), with inactive but functional SEM (self-elimination mechanism) element
#'  \item U: Low-cost resistant allele, non-targetable by GD or SEM (R1 in other literature)
#'  \item R: High-cost resistant allele, non-targetable by GD or SEM (R2 in other literature)
#'  \item V: "wild-type", product of SEM, non-targetable by GD or SEM
#'  \item H: Active GD with active SEM
#'  \item S: Active GD with non-functional SEM
#'  \item Y: Y-chromosome, non-targetable
#' }
#'
#' "V" alleles are simply a minor allele with the same protein sequence as the
#' major allele, "W". There is the possibility for allelic conversion of the "V"
#' allele into the "W" allele by mechanisms such as MMR.
#'
#' This drive has female specific GD and SEM parameters, as well as \strong{correct}
#' maternal deposition. There are no dosage effects
#' modeled (i.e., having two GD alleles increasing or decreasing the GD rates).
#'
#'
#' @param pF Rate of cleavage during GD process in females
#' @param qF Rate of HDR during GD process in females
#' @param rF Rate of in-frame resistance generation during GD process in females
#'
#' @param aF Rate of cleavage during SEM process in females
#' @param bF Rate of SSA during SEM process in females
#' @param cF Rate of "V" allele formation from SSA during SEM process in females
#'
#' @param mmrF Rate of MMR in females, driving allelic conversion of "V" into "W"
#' @param mmrM Rate of MMR in males, driving allelic conversion of "V" into "W"
#'
#' @param pDep Rate of cleavage during maternal deposition
#' @param qDep Rate of HDR during maternal deposition
#' @param rDep Rate of in-frame resistance generation during maternal deposition
#'
#' @param eta Genotype-specific mating fitness
#' @param phi Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF Genotype-specific female pupatory success
#' @param xiM Genotype-specific male pupatory success
#' @param s Genotype-specific fractional reduction(increase) in fertility
#'
#'
#' @return Named list containing the inheritance cube, transition matrix, genotypes,
#' wild-type allele, and all genotype-specific parameters.
#' @export
cubeSEMX <- function(pF=1, qF=1, rF=0,
                     aF=1, bF=1, cF=1,
                     mmrF=0, mmrM=0,
                     pDep=0, qDep=0, rDep=0,
                     eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  ## safety checks
  inputVec <- c(pF,qF,rF, aF,bF,cF, mmrF,mmrM, pDep,qDep,rDep)
  if(any(inputVec>1) || any(inputVec<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }


  # # Testing Probs
  # testVec <- runif(n = 19, min = 0, max = 1)
  #
  # mmrF <- testVec[1]; mmrM <- testVec[2]
  # pF <- testVec[3]; qF <- testVec[4]; rF <- testVec[5];
  # aF <- testVec[6]; bF <- testVec[7]; cF <- testVec[8]; dF <- testVec[9];
  # pDep <- testVec[17]; qDep <- testVec[18]; rDep <- testVec[19];


  #############################################################################
  ## generate genotypes
  #############################################################################

  # # List of possible alleles
  # #  This time, I'll split male/female, so that we can tack them together at the
  # #  end. This will provide a nice organization in the resulting cube.
  # alleles <- c('X', 'G', 'U', 'R', 'V', 'H', 'S')
  # # Generate genotypes
  # # This generates all combinations of alleles
  # holdF <- as.matrix(x = expand.grid(alleles, alleles, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  # holdM <- as.matrix(x = expand.grid(alleles, 'Y', KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  # # Sort both alleles to remove uniques
  # #  The alleles aren't unique, and can occur in any order, so we sort them
  # #  alphabetically in order to remove duplicate genotypes
  # sortedAllelesF <- apply(X = holdF, MARGIN = 1, FUN = sort)
  # sortedAllelesM <- apply(X = holdM, MARGIN = 1, FUN = sort)
  # # Paste together and get unique ones only
  # genoF <- unique(do.call(what = paste0, args = list(sortedAllelesF[1, ], sortedAllelesF[2, ])))
  # genoM <- unique(do.call(what = paste0, args = list(sortedAllelesM[1, ], sortedAllelesM[2, ])))

  genoF <- c('XX','GX','UX','RX','VX','HX','SX',
             'GG','GU','GR','GV','GH','GS','UU',
             'RU','UV','HU','SU','RR','RV','HR',
             'RS','VV','HV','SV','HH','HS','SS')
  genoM <- c('XY','GY','UY','RY','VY','HY','SY')
  genotypes <- c(genoF, genoM)


  #############################################################################
  ## setup probability lists
  #############################################################################
  # This assumes that each process goes to completion before another process begins
  #  i.e., if there's GD, then only GD happens
  #  but, if there's SEM, then GD happens first, and then SEM is performed on the
  #   results of the GD process.
  # IF THESE DON'T OCCUR IN THIS ORDER, THE PROBABILITIES BELOW WILL BE INCCORECT.
  #
  # For males, it is the same assumptions for GD then SEM.
  # After that, the possible male alleles are matched with possible female alleles,
  #  and then maternal deposition is applied. This allows for proper HDR during
  #  deposition.
  #
  # Remember, "H" is never inherited, it's turned on, so "H" becomes "G" always

  ##########
  # Female
  ##########
  mendF <- list('X' = c('X'=1),
                'U' = c('U'=1),
                'R' = c('R'=1),
                'V' = c('X'=mmrF, 'V'=1-mmrF))

  # assume that G, H, and S are equally competent at homing
  gdF <- list('X' = list('G' = c('X'=1-pF,
                                 'G'=pF*qF,
                                 'U'=pF*(1-qF)*rF,
                                 'R'=pF*(1-qF)*(1-rF)),
                         'H' = c('X'=1-pF,
                                 'G'=pF*qF,
                                 'U'=pF*(1-qF)*rF,
                                 'R'=pF*(1-qF)*(1-rF)),
                         'S' = c('X'=1-pF,
                                 'S'=pF*qF,
                                 'U'=pF*(1-qF)*rF,
                                 'R'=pF*(1-qF)*(1-rF)) ),
              'G' = c('G'=1),
              'U' = c('U'=1),
              'R' = c('R'=1),
              'V' = c('X'=mmrF, 'V'=1-mmrF),
              'H' = c('G'=1),
              'S' = c('S'=1))

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "X" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  #  "H" should have been turned into "G" by here, which is fine since the
  #      outcomes are the same
  semF <- list('X' = c('X'=1),
               'G' = c('G'=1-aF,
                       'V'=aF*bF*cF,
                       'X'=aF*bF*(1-cF),
                       'S'=aF*(1-bF) ),
               'U' = c('U'=1),
               'R' = c('R'=1),
               'V' = c('V'=1),
               'H' = c('G'=1-aF,
                       'V'=aF*bF*cF,
                       'X'=aF*bF*(1-cF),
                       'S'=aF*(1-bF) ),
               'S' = c('S'=1))

  ##########
  # Male
  ##########
  mendM <- list('X' = c('X'=1),
                'G' = c('G'=1),
                'U' = c('U'=1),
                'R' = c('R'=1),
                'V' = c('X'=mmrM, 'V'=1-mmrM),
                'H' = c('G'=1),
                'S' = c('S'=1),
                'Y' = c('Y'=1))

  # only doing GD deposition - can assume that SEM is controlled and not under
  #  germline expression.
  #  GD deposition only impacts the 'X' allele, so it's the only one we need to
  #  worry about
  # Since we don't know what happened in the female before we got here, we can
  #  have nearly any allele to repair against. The only one we shouldn't see is
  #  "H", so I'll leave it out, and if we do see it, figure out why.
  # Since this is HDR dependent, I give the "V" one an opportunity for increased
  #  MMR, thereby double converting new "V" alleles on into "X" alleles
  dep <- list('X' = c('X'=1-pDep + pDep*qDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'G' = c('X'=1-pDep,
                      'G'=pDep*qDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'U' = c('X'=1-pDep,
                      'U'=pDep*qDep + pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'R' = c('X'=1-pDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*qDep + pDep*(1-qDep)*(1-rDep)),
              'V' = c('X'=1-pDep + pDep*qDep*mmrM,
                      'V'=pDep*qDep*(1-mmrM),
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'S' = c('X'=1-pDep,
                      'S'=pDep*qDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)) )


  #############################################################################
  ## Fill transition matrix
  #############################################################################

  # Use this many times down below
  numGenF <- length(genoF)
  numGenM <- length(genoM)
  numGen <- numGenF + numGenM

  # Create transition matrix to fill
  tMatrix <- array(data = 0,
                   dim = c(numGen,numGen,numGen),
                   dimnames = list(genotypes,genotypes,genotypes))


  #############################################################################
  ## Loop over all matings, female outer loop
  #############################################################################
  for(fi in 1:numGenF){
    # Female loop
    # Split genotype into alleles
    fSplit <- strsplit(x = genotypes[fi], split = "", useBytes = TRUE)[[1]]

    # Score them
    semScoreF <- any("H" == fSplit)
    gdScoreF <- semScoreF || any("G" == fSplit) || any("S" == fSplit)


    ##########
    # Female Alleles
    ##########
    # First step
    if(!gdScoreF){
      # Mendelian
      fAllele <- c(mendF[[ fSplit[1] ]], mendF[[ fSplit[2] ]])

    } else {
      # GD

      # check what the second allele is
      #  If "W", it has a different shape, which depends on which GD allele is
      #  present. Since "W" is only in the second place, we can check it easily.
      if(fSplit[[2]]=='X'){
        # There is a "W" in the second spot
        fAllele <- c(gdF[[ fSplit[1] ]], gdF[['X']][[ fSplit[1] ]])
      } else {
        # There is not a "W" in the second spot
        fAllele <- c(gdF[[ fSplit[1] ]], gdF[[ fSplit[2] ]])
      }

    } # end first step, Mendelian vs GD

    # Second Step
    if(semScoreF){
      # SEM
      # Here, we take each allele from above, and we target it with the SEM construct
      #  This is done by multiplying the breakdown of GD results by possible SEM results.
      # Finally, unlist so it is the same shape as the previous two.
      # "Map()" is faster, but I can't figure out the output naming.
      fAllele <- unlist(x = mapply(FUN = "*", fAllele, semF[names(fAllele)],
                                   SIMPLIFY = FALSE, USE.NAMES = FALSE),
                        recursive = TRUE, use.names = TRUE)
    }

    # Aggregate duplicate alleles
    #  This is mostly useful for SEM, and I could have included it in that step.
    #  However, it will also be often useful in the GD step (but not always)
    #  and is useful in one case in the Mendelian step
    fAlleleReduc <- vapply(X = unique(names(fAllele)),
                           FUN = function(x){sum(fAllele[names(fAllele)==x])},
                           FUN.VALUE = numeric(length = 1L))


    ###########################################################################
    ## Male loop. This is the inner loop
    ###########################################################################
    for(mi in (1:numGenM + numGenF)){
      # Male loop
      # Split genotype into alleles
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]

      # Score them
      depScoreM <- any("X" == mSplit)


      ##########
      # Male Alleles
      ##########
      # Always Mendelian, since there's nothing for the drive to target in males
      mAlleleReduc <- c(mendM[[ mSplit[1] ]], mendM[[ mSplit[2] ]])

      # # I don't think this is necessary in this cube
      # # Aggregate duplicate alleles
      # #  This is mostly useful for SEM, and I could have included it in that step.
      # #  However, it will also be often useful in the GD step (but not always)
      # #  and is useful in one case in the Mendelian step
      # mAlleleReduc <- vapply(X = unique(names(mAllele)),
      #                        FUN = function(x){sum(mAllele[names(mAllele)==x])},
      #                        FUN.VALUE = numeric(length = 1L))


      #########################################################################
      ## Get combinations of male/female alleles
      ##  Perform maternal deposition if relevant
      ## Put results into tMatrix
      #########################################################################
      # The maternal deposition check and loop is not very performant, and has
      #  significant code duplication. I can't think of a better way to do it without
      #  setting up an awkward list for deposition alleles. So, we duplicate the
      #  probability and allele sorting/pasting code into 3 places - no deposition,
      #  deposition with no allele impacted, and deposition against "W" alleles.
      #  Additionally, we grow a vector in the deposition section, which isn't great
      #  but I don't have a better idea (and it's small, so hopefully not terrible).

      ##########
      # Maternal deposition check
      ##########
      #  if there are any 'W' alleles in the male germline, and any GD alleles
      #  from the females
      if(gdScoreF && depScoreM){
        # There is meaningful deposition

        # return object
        genoVec <- numeric(length = 0L)
        genoNames <- character(length = 0L)

        # things I reference a lot
        mARNames <- names(mAlleleReduc)
        fARNames <- names(fAlleleReduc)

        # loop over male alleles
        for(eleM in 1:length(mAlleleReduc)){
          # check for allele type
          #  only "W" alleles are impacted by deposition, so only they have
          #  complicated repair structures
          if(mARNames[eleM] == 'X'){
            # "W" allele

            # Repair during deposition depends on the maternal allele inherited,
            #  so now we loop over all female alleles, and use that to define
            #  the repair process.
            # The we process everything and append to geno/probs vectors
            for(eleF in 1:length(fAlleleReduc)){
              # Get repair distribution
              repVec <- dep[[ fARNames[eleF] ]]

              ##########
              # Genotype probs
              ##########
              #  These come from the male allele probs, the female allele probs,
              #  and the distribution of repair products
              genoVec <- c(genoVec, mAlleleReduc[eleM] * fAlleleReduc[eleF] * repVec)

              ##########
              # Genotypes
              ##########
              #  Get all combinations of alleles
              oneAGenos <- expand.grid(names(repVec), fARNames[eleF],
                                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
              #  Sort alleles
              oneAGenos <- apply(X = oneAGenos, MARGIN = 1, FUN = sort.int)
              #  paste and append genotypes
              genoNames <- c(genoNames, file.path(oneAGenos[1,], oneAGenos[2,], fsep = ""))

            } # end female allele loop

          } else {
            # Not a "W" allele
            # These alleles are NOT impacted by deposition, so their outcome
            #  is simply the probabiliy of the male allele, multiplied by the
            #  probability of each female allele

            ##########
            # Genotype probs
            ##########
            genoVec <- c(genoVec, mAlleleReduc[eleM] * fAlleleReduc)

            ##########
            # Genotypes
            ##########
            #  Get all combinations of alleles
            oneAGenos <- expand.grid(mARNames[eleM], fARNames, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
            #  Sort alleles
            oneAGenos <- apply(X = oneAGenos, MARGIN = 1, FUN = sort.int)
            #  paste and append genotypes
            genoNames <- c(genoNames, file.path(oneAGenos[1,], oneAGenos[2,], fsep = ""))

          } # end male allele check

        } # end male allele loop


      } else {
        # There is not meaningful deposition
        # Alleles recombine independently
        holdProbs <- expand.grid(fAlleleReduc, mAlleleReduc,
                                 KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele <- expand.grid(names(fAlleleReduc), names(mAlleleReduc),
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # sort alleles into order
        holdAllele <- apply(X = holdAllele, MARGIN = 1, FUN = sort.int)

        ##########
        # Genotype Probs
        ##########
        genoVec <- holdProbs[ ,1]*holdProbs[ ,2]

        ##########
        # Genotypes
        ##########
        genoNames <- file.path(holdAllele[1,], holdAllele[2,], fsep = "")

      } # end maternal deposition and allele -> genotypes calculation


      ##########
      # Finish Genotype Distribution
      ##########
      # Aggregate duplicate genotypes
      genoAgg <- vapply(X = unique(genoNames),
                        FUN = function(x){sum(genoVec[genoNames==x])},
                        FUN.VALUE = numeric(length = 1L))

      #normalize
      genoAgg <- genoAgg/sum(genoAgg)

      #set values in tMatrix
      tMatrix[fi,mi, names(genoAgg) ] <- genoAgg


    }# end male loop
  }# end female loop


  # # Test stuff
  # for(mi in 1:numGen){
  #   for(fi in 1:numGen){
  #     if(sum(tMatrix[fi,mi, ])-1 > 0.0008){
  #       cat("Male:", genotypes[mi], ", female:", genotypes[fi])
  #     }
  #   }
  # }


  ## Protect from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0

  ## Initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1.0, dim = c(numGen,numGen,numGen),
                         dimnames = list(genotypes, genotypes, genotypes))

  ## Genotype-specific modifiers
  modifiers = cubeModifiers(gtype = genotypes, eta = eta, phi = phi,
                            omega = omega, xiF = xiF, xiM = xiM, s = s)


  ## Put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = genotypes,
    genotypesN = numGen,
    wildType = c('XX','XY'),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = 'GY'
  ))

}
