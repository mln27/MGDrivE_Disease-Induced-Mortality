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


  # Testing Probs
  testVec <- runif(n = 16, min = 0, max = 1)

  mmrF <- testVec[1]; mmrM <- testVec[2]
  pF <- testVec[3]; qF <- testVec[4]; rF <- testVec[5];
  aF <- testVec[6]; bF <- testVec[7]; cF <- testVec[8]; dF <- testVec[9];
  pM <- testVec[10]; qM <- testVec[11]; rM <- testVec[12];
  aM <- testVec[13]; bM <- testVec[14]; cM <- testVec[15]; dM <- testVec[16];









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
  #  GD deposition only impacts the 'W' allele, so it's the only one we need to
  #  worry about




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




  c('W', 'G', 'U', 'R', 'V', 'H', 'S')










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
  ## Loop over all matings, female outer loop
  #############################################################################
  for(fi in 1:numGen){
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
      if(fSplit[[2]]=='W'){
        # There is a "W" in the second spot
        fAllele <- c(gdF[[ fSplit[1] ]], gdF[['W']][[ fSplit[1] ]])
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
                                   SIMPLIFY = TRUE, USE.NAMES = FALSE),
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
    for(mi in 1:numGen){
      # Male loop
      # Split genotype into alleles
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]

      # Score them
      semScoreM <- any("H" == mSplit)
      gdScoreM <- semScoreM || any("G" == mSplit) || any("S" == mSplit)
      depScoreM <- any("W" == mSplit)


      ##########
      # Male Alleles
      ##########
      # First step
      if(!gdScoreM){
        # Mendelian
        mAllele <- c(mendM[[ mSplit[1] ]], mendM[[ mSplit[2] ]])

      } else {
        # GD

        # check what the second allele is
        #  If "W", it has a different shape, which depends on which GD allele is
        #  present. Since "W" is only in the second place, we can check it easily.
        if(mSplit[[2]]=='W'){
          # There is a "W" in the second spot
          mAllele <- c(gdM[[ mSplit[1] ]], gdM[['W']][[ mSplit[1] ]])
        } else {
          # There is not a "W" in the second spot
          mAllele <- c(gdM[[ mSplit[1] ]], gdM[[ mSplit[2] ]])
        }

      } # end first step, Mendelian vs GD

      # Second Step
      if(semScoreM){
        # SEM
        # Here, we take each allele from above, and we target it with the SEM construct
        #  This is done by multiplying the breakdown of GD results by possible SEM results.
        # Finally, unlist so it is the same shape as the previous two.
        # "Map()" is faster, but I can't figure out the output naming.
        mAllele <- unlist(x = mapply(FUN = "*", mAllele, semM[names(mAllele)],
                                     SIMPLIFY = TRUE, USE.NAMES = FALSE),
                          recursive = TRUE, use.names = TRUE)
      }

      # Aggregate duplicate alleles
      #  This is mostly useful for SEM, and I could have included it in that step.
      #  However, it will also be often useful in the GD step (but not always)
      #  and is useful in one case in the Mendelian step
      mAlleleReduc <- vapply(X = unique(names(mAllele)),
                             FUN = function(x){sum(mAllele[names(mAllele)==x])},
                             FUN.VALUE = numeric(length = 1L))


      #########################################################################
      ## Get combinations of male/female alleles
      ##  Perform maternal deposition if relevant
      ## Put results into tMatrix
      #########################################################################
      # Maternal deposition check
      #  if there are any 'W' alleles in the male germline, and any GD alleles
      #  from the females
      if(gdScoreF && depScoreM){
        # There is meaningful deposition

        holdList <- vector(mode = "list", length = length(mAlleleReduc))

        for(elem in 1:length(mAlleleReduc)){
          holdList[[elem]] <-
        }












      } else {
        # There is not meaningful deposition
        # Alleles recombine indepdendently
        holdProbs <- expand.grid(fAlleleReduc, mAlleleReduc,
                                 KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele <- expand.grid(names(fAlleleReduc), names(mAlleleReduc),
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)










      }








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
