###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE2: Mosquito Gene Drive Explorer V2 with Disease-Induced Mortality
#   ReMEDE (REpeat Mediated Excision of a Drive Element)
#   Jared Bennett
#   jared_bennett@berkeley.edu
#   March 2022
#
###############################################################################

#' Inheritance Cube: ReMEDE (REpeat Mediated Excision of a Drive Element)
#'
#' The ReMEDE system, put forth by [Chennuri and Myles](GETLINKWHENEXISTS!!!),
#' is a realized application of [biodegradable gene drives](https://doi.org/10.1098/rstb.2019.0804)
#' with application of [Small-Molecule Control](https://doi.org/10.1016/j.celrep.2020.107841)
#' for the removal of drive elements through induceable SSA. This construct has a
#' single, autosomal target site that carries gRNA, CAS9, an induceable endonuclease,
#' and is flanked by direct repeats. There are 7 possible alleles:
#'  * W: Wild-type
#'  * G: Active (GD) gene drive, with inactive but functional SEM (self-elimination mechanism) element
#'  * U: Low-cost resistant allele, non-targetable by GD or SEM (R1 in other literature)
#'  * R: High-cost resistant allele, non-targetable by GD or SEM (R2 in other literature)
#'  * V: "wild-type", product of SEM, non-targetable by GD or SEM
#'  * H: Active GD with active SEM
#'  * S: Active GD with non-functional SEM
#'
#' "V" alleles are simply a minor allele with the same protein sequence as the
#' major allele, "W". There is the possibility for allelic conversion of the "V"
#' allele into the "W" allele by mechanisms such as MMR.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' This drive has male and female specific GD and SEM parameters, as well as
#' PROPER FEMALE DEPOSITION???????? There are no dosage effets modeled (i.e., having
#' two GD alleles increasing or decreasing the GD rates).
#'
#' @param mmrF Rate of MMR in females, driving allelic conversion of "V" into "W"
#' @param mmrM Rate of MMR in males, driving allelic conversion of "V" into "W"
#'
#' @param pF Rate of cleavage during GD process in females
#' @param qF Rate of HDR during GD process in females
#' @param rF Rate of in-frame resistance generation during GD process in females
#'
#' @param aF Rate of cleavage during SEM process in females
#' @param bF Rate of SSA during SEM process in females
#' @param cF Rate of "V" allele formation from SSA during SEM process in females
#' @param dF Rate of "S" vs "R" allele formation from NHEJ durin SEM process in females
#'
#'
#' @param pM Rate of cleavage during GD process in males
#' @param qM Rate of HDR during GD process in males
#' @param rM Rate of in-frame resistance generation during GD process in males
#'
#' @param aM Rate of cleavage during SEM process in males
#' @param bM Rate of SSA during SEM process in males
#' @param cM Rate of "V" allele formation from SSA during SEM process in males
#' @param dM Rate of "S" vs "R" allele formation from NHEJ durin SEM process in males
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @param cM1 Maternally inherited Cas9 cutting rate at locus 1
#' @param cM2 Maternally inherited Cas9 cutting rate at locus 2
#' @param cP1 Paternally inherited Cas9 cutting rate at locus 1
#' @param cP2 Paternally inherited Cas9 cutting rate at locus 2
#' @param hM1 Maternally inherited Cas9 homing efficiency at locus 1
#' @param hM2 Maternally inherited Cas9 homing efficiency at locus 2
#' @param hP1 Paternally inherited Cas9 homing efficiency at locus 1
#' @param hP2 Paternally inherited Cas9 homing efficiency at locus 2
#' @param rM1 Maternally inherited Cas9 resistance efficiency at locus 1
#' @param rM2 Maternally inherited Cas9 resistance efficiency at locus 2
#' @param rP1 Paternally inherited Cas9 resistance efficiency at locus 1
#' @param rP2 Paternally inherited Cas9 resistance efficiency at locus 2
#' @param crM Maternal crossover rate, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crP Paternal crossover rate, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#'
#'
#'
#'
#'
#' @param eta Genotype-specific mating fitness
#' @param phi Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF Genotype-specific female pupatory success
#' @param xiM Genotype-specific male pupatory success
#' @param s Genotype-specific fractional reduction(increase) in fertility
#'
#'
#'
#'
#'
#'
#'
#' @return Named list containing the inheritance cube, transition matrix, genotypes, wild-type allele,
#' and all genotype-specific parameters.
#' @export
cubeSEM <- function(cM1=0, cM2=0, cP1=0, cP2=0,
                     hM1=0, hM2=0, hP1=0, hP2=0,
                     rM1=0, rM2=0, rP1=0, rP2=0, crM=0, crP=0,
                     eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  ## safety checks
  if(any(c(cM1,cM2,cP1,cP2,hM1,hM2,hP1,hP2,rM1,rM2,rP1,rP2,crM,crP)>1) || any(c(cM1,cM2,cP1,cP2,hM1,hM2,hP1,hP2,rM1,rM2,rP1,rP2,crM,crP)<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }












  #############################################################################
  ## generate genotypes
  #############################################################################

  # # List of possible alleles
  # alleles <- c('W', 'G', 'U', 'R', 'V', 'H', 'S')
  # # Generate genotypes
  # # This generates all combinations of alleles
  # hold <- as.matrix(x = expand.grid(alleles, alleles, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  # # Sort both alleles to remove uniques
  # #  The alleles aren't unique, and can occur in any order, so we sort them
  # #  alphabetically in order to remove duplicate genotypes
  # sortedAlleles <- apply(X = hold, MARGIN = 1, FUN = sort)
  # # Paste together and get unique ones only
  # genotypes <- unique(do.call(what = paste0, args = list(sortedAlleles[1, ], sortedAlleles[2, ])))

  genotypes <- c('WW','GW','UW','RW','VW','HW','SW',
                 'GG','GU','GR','GV','GH','GS','UU',
                 'RU','UV','HU','SU','RR','RV','HR',
                 'RS','VV','HV','SV','HH','HS','SS')


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
  mendF <- list('W' = c('W'=1),
                'U' = c('U'=1),
                'R' = c('R'=1),
                'V' = c('W'=mmrF, 'V'=1-mmrF))

  # assume that G, H, and S are equally competent at homing
  gdF <- list('W' = list('G' = c('W'=1-pF,
                                 'G'=pF*qF,
                                 'U'=pF*(1-qF)*rF,
                                 'R'=pF*(1-qF)*(1-rF)),
                         'H' = c('W'=1-pF,
                                 'G'=pF*qF,
                                 'U'=pF*(1-qF)*rF,
                                 'R'=pF*(1-qF)*(1-rF)),
                         'S' = c('W'=1-pF,
                                 'S'=pF*qF,
                                 'U'=pF*(1-qF)*rF,
                                 'R'=pF*(1-qF)*(1-rF)) ),
              'G' = c('G'=1),
              'U' = c('U'=1),
              'R' = c('R'=1),
              'V' = c('W'=mmrF, 'V'=1-mmrF),
              'H' = c('G'=1),
              'S' = c('S'=1))

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "W" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  #  "H" should have been turned into "G" by here, which is fine since the
  #      outcomes are the same
  semF <- list('W' = c('W'=1),
               'G' = c('G'=1-aF,
                       'V'=aF*bF*cF,
                       'W'=aF*bF*(1-cF),
                       'S'=aF*(1-bF)*dF,
                       'R'=aF*(1-bF)*(1-dF) ),
               'U' = c('U'=1),
               'R' = c('R'=1),
               'V' = c('V'=1),
               'H' = c('G'=1-aF,
                       'V'=aF*bF*cF,
                       'W'=aF*bF*(1-cF),
                       'S'=aF*(1-bF)*dF,
                       'R'=aF*(1-bF)*(1-dF) ),
               'S' = c('S'=1))

  ##########
  # Male
  ##########
  mendM <- list('W' = c('W'=1),
                'U' = c('U'=1),
                'R' = c('R'=1),
                'V' = c('W'=mmrM, 'V'=1-mmrM))

  # assume that G, H, and S are equally competent at homing
  gdM <- list('W' = list('G' = c('W'=1-pM,
                                 'G'=pM*qM,
                                 'U'=pM*(1-qM)*rM,
                                 'R'=pM*(1-qM)*(1-rM)),
                         'H' = c('W'=1-pM,
                                 'G'=pM*qM,
                                 'U'=pM*(1-qM)*rM,
                                 'R'=pM*(1-qM)*(1-rM)),
                         'S' = c('W'=1-pM,
                                 'S'=pM*qM,
                                 'U'=pM*(1-qM)*rM,
                                 'R'=pM*(1-qM)*(1-rM)) ),
              'G' = c('G'=1),
              'U' = c('U'=1),
              'R' = c('R'=1),
              'V' = c('W'=mmrM, 'V'=1-mmrM),
              'H' = c('G'=1),
              'S' = c('S'=1))

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "W" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  #  "H" should have been turned into "G" by here, which is fine since the
  #      outcomes are the same
  semM <- list('W' = c('W'=1),
               'G' = c('G'=1-aM,
                       'V'=aM*bM*cM,
                       'W'=aM*bM*(1-cM),
                       'S'=aM*(1-bM)*dM,
                       'R'=aM*(1-bM)*(1-dM) ),
               'U' = c('U'=1),
               'R' = c('R'=1),
               'V' = c('V'=1),
               'H' = c('G'=1-aM,
                       'V'=aM*bM*cM,
                       'W'=aM*bM*(1-cM),
                       'S'=aM*(1-bM)*dM,
                       'R'=aM*(1-bM)*(1-dM) ),
               'S' = c('S'=1))

  # only doing GD deposition - can assume that SEM is controlled and not under
  #  germline expression.
  # There shouldn't be any 'H' by the time we hit this, but in case there is,
  #  inherit as 'G'
  # Don't let 'V' undergo a second round of MMR
  dep <- list('W' = list('W', 'G', 'U', 'R', 'V', 'H', 'S'),





              'G' = c('G'=1),
              'U' = c('U'=1),
              'R' = c('R'=1),
              'V' = c('V'=1),
              'H' = c('G'=1),
              'S' = c('S'=1))













  #############################################################################
  ## Fill transition matrix
  #############################################################################

  # Use this many times down below
  numGen <- length(genotypes)

  # Create transition matrix to fill
  tMatrix <- array(data = 0,
                   dim = c(numGen,numGen,numGen),
                   dimnames = list(genotypes,genotypes,genotypes))


  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for(fi in 1:numGen){
    # Female loop
    # Split genotype into alleles
    fSplit <- strsplit(x = genotypes[fi], split = "", useBytes = TRUE)[[1]]

    # Score them
    #  Because R is fun, we can check and store the value of the H check, while
    #  also using it in the expression for G and S
    gdScoreF <- any("G" == fSplit) || (semScoreF <- any("H" == fSplit)) || any("S" == fSplit)














    #setup offspring allele lists
    fPHold <- rep(x = list(list()), times = numAlleles) #vector(mode = "list", length = numAlleles)
    fAllele <- character(length = 0L)
    fProbs <- numeric(length = 0L)


    # check if there is a Cas9 from either parent
    if(fScore[1] && fScore[2]){
      # There is a Cas9 from one parent or the other
      # There are gRNAs
      # Thus, there is homing

      # loop over alleles, must keep target sites linked
      for(allele in 1:numAlleles){
        # we make the assumption that female homing is more important than male
        # ie, if we have M and P at locus one on both alleles, the M takes precedence
        #  in terms of homing.

        # check if maternal homing
        if(any(fAlleleMat[ ,1]=="M")){
          # at least one of the Cas9 proteins is from the mother
          fPHold[[allele]][[1]] <- fLocus1$maternal[[ fAlleleMat[allele,1] ]]
          fPHold[[allele]][[2]] <- Locus2$maternal[[ fAlleleMat[allele,2] ]]
        } else {
          # the Cas9 protein is from the father
          fPHold[[allele]][[1]] <- fLocus1$paternal[[ fAlleleMat[allele,1] ]]
          fPHold[[allele]][[2]] <- Locus2$paternal[[ fAlleleMat[allele,2] ]]
        } # end M/P homing

      } # end loop over alleles
    } else {
      # either no Cas9 or no gRNAs
      # all inheritance is Mendelian

      # loop over alleles, must keep target sites linked
      for(allele in 1:numAlleles){
        # set target site 1, then target site 2
        fPHold[[allele]][[1]] <- fLocus1$mendelian[[ fAlleleMat[allele,1] ]]
        fPHold[[allele]][[2]] <- Locus2$mendelian[[ fAlleleMat[allele,2] ]]

      } # end loop over alleles
    } # end female scoring


    # perform cross-overs
    hold1 <- fPHold[[1]][[1]] # need

    fPHold[[1]][[1]] <- c((1-crM)*hold1, crM*fPHold[[2]][[1]])
    fPHold[[2]][[1]] <- c((1-crM)*fPHold[[2]][[1]], crM*hold1)

    # all combinations of female alleles.
    for(allele in 1:numAlleles){
      # make combinations of the allele, then store those combinations to mix
      #  with the next allele
      # expand combinations
      holdProbs <- expand.grid(fPHold[[allele]],KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(names(fPHold[[allele]][[1]]), names(fPHold[[allele]][[2]]),
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # contract (paste or multiply) combinations and store
      fProbs <- c(fProbs, holdProbs[,1]*holdProbs[,2])
      fAllele <- c(fAllele, file.path(holdAllele[,1], holdAllele[,2], fsep = ""))
    }

    # remove zeros
    fAllele <- fAllele[fProbs!=0]
    fProbs <- fProbs[fProbs!=0]


    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    for (mi in 1:numGen){
      # isn't symmetric

      #do male stuff here
      #split male genotype
      #This splits all characters.
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]
      #Matrix of alleles and target sites
      #  Each row is an allele, there are 2 alleles
      #  Each column is the target site on that allele, there are 2 target sites
      #  this is because target sites are linked on an allele
      mAlleleMat <- matrix(data = mSplit, nrow = numAlleles, ncol = numAlleles, byrow = TRUE)
      #Score them
      mScore[1] <- any("P" == mAlleleMat) || any("M" == mAlleleMat)
      mScore[2] <- any("G" == mAlleleMat)

      #setup offspring allele lists
      mPHold <- rep(x = list(list()), times = numAlleles) #vector(mode = "list", length = numAlleles)
      mAllele <- character(length = 0L)
      mProbs <- numeric(length = 0L)


      # check if there is a Cas9 from either parent
      if(mScore[1] && mScore[2]){
        # There is a Cas9 from one parent or the other
        # There are gRNAs
        # Thus, there is homing

        # loop over alleles, must keep target sites linked
        for(allele in 1:numAlleles){
          # we make the assumption that female homing is more important than male
          # ie, if we have M and P at locus one on both alleles, the M takes precedence
          #  in terms of homing.

          # check if maternal homing
          if(any(mAlleleMat[ ,1]=="M")){
            # at least one of the Cas9 proteins is from the mother
            mPHold[[allele]][[1]] <- mLocus1$maternal[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- Locus2$maternal[[ mAlleleMat[allele,2] ]]
          } else {
            # the Cas9 protein is from the father
            mPHold[[allele]][[1]] <- mLocus1$paternal[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- Locus2$paternal[[ mAlleleMat[allele,2] ]]
          } # end M/P homing

        } # end loop over alleles
      } else {
        # either no Cas9 or no gRNAs
        # all inheritance is Mendelian

        # loop over alleles, must keep target sites linked
        for(allele in 1:numAlleles){
          # set target site 1, then target site 2
          mPHold[[allele]][[1]] <- mLocus1$mendelian[[ mAlleleMat[allele,1] ]]
          mPHold[[allele]][[2]] <- Locus2$mendelian[[ mAlleleMat[allele,2] ]]

        } # end loop over alleles
      } # end male scoring


      # perform cross-overs
      hold1 <- mPHold[[1]][[1]] # need

      mPHold[[1]][[1]] <- c((1-crP)*hold1, crP*mPHold[[2]][[1]])
      mPHold[[2]][[1]] <- c((1-crP)*mPHold[[2]][[1]], crP*hold1)

      # all combinations of female alleles.
      for(allele in 1:numAlleles){
        # make combinations of the allele, then store those combinations to mix
        #  with the next allele
        # expand combinations
        holdProbs <- expand.grid(mPHold[[allele]],KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele <- expand.grid(names(mPHold[[allele]][[1]]), names(mPHold[[allele]][[2]]),
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # contract (paste or multiply) combinations and store
        mProbs <- c(mProbs, holdProbs[,1]*holdProbs[,2])
        mAllele <- c(mAllele, file.path(holdAllele[,1], holdAllele[,2], fsep = ""))
      }

      # remove zeros for test
      mAllele <- mAllele[mProbs!=0]
      mProbs <- mProbs[mProbs!=0]


      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################

      # male and female alleles/probs are allready combined by target sites, and
      #  we assume alleles segregate indepenently, so we just have to get combinations
      #  of male vs female for the offspring
      holdProbs <- expand.grid(fProbs, mProbs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(fAllele, mAllele, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # sort alleles into order
      holdAllele <- apply(X = holdAllele, MARGIN = 1, FUN = sort.int)

      # contract (paste or multiply) combinations
      holdProbs <- holdProbs[ ,1]*holdProbs[ ,2]
      holdAllele <- file.path(holdAllele[1,], holdAllele[2,], fsep = "")

      #aggregate duplicate genotypes
      reducedP <- vapply(X = unique(holdAllele), FUN = function(x){
        sum(holdProbs[holdAllele==x])},
        FUN.VALUE = numeric(length = 1L))

      #normalize
      reducedP <- reducedP/sum(reducedP)

      #set values in tMatrix
      tMatrix[fi,mi, names(reducedP) ] <- reducedP

    }# end male loop
  }# end female loop

  ## set the other half of the matrix
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors

  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(numGen,numGen,numGen),
                         dimnames = list(genotypes, genotypes, genotypes))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(genotypes, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = genotypes,
    genotypesN = numGen,
    wildType = "WWWW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "MGPG"
  ))

}
