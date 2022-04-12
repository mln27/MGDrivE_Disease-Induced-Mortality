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

#' Inheritance Cube: ReMEDE (REpeat Mediated Excision of a Drive Element) in Trans
#'
#' This is the second ReMEDE system put forth by \href{GETLINKWHENEXISTS!!!}{Chennuri and Myles}.
#' It's another realization of \href{https://doi.org/10.1098/rstb.2019.0804}{biodegradable gene drives},
#' where the pieces exist as indeptendently-recombining constructs that only act
#' when in the presence of each other.
#'
#' This construct consists of two autosomal target sites. The first site carries
#' gRNA and Cas9 for homing and removal of the SEM site. The second site consists
#' of an endonuclease targeting the Cas9 site. Both sites contain direct repeats
#' to induce SSA upon chromosome damage. There are 6 possible alleles at the first site:
#' \itemize{
#'  \item W: Wild-type allele
#'  \item G: Active GD (gene drive), with targetable SEM (self-elimination mechanism) site
#'  \item U: Low-cost resistant allele, non-targetable by GD or SEM (R1 in other literature)
#'  \item R: High-cost resistant allele, non-targetable by GD or SEM (R2 in other literature)
#'  \item V: "wild-type", product of SEM, non-targetable by GD or SEM
#'  \item S: Active GD with non-targetable SEM site
#' }
#'
#' There are 4 possible alleles at the second site:
#' \itemize{
#'  \item W: Wild-type allele
#'  \item H: Active SEM element
#'  \item R: High-cost resistant allele, non-targetable by GD (R2 in other literature)
#'  \item E: "wild-type", produce of GD, non-targetbale by GD
#' }
#'
#' This provides a total genotype count of 210. Genotypes are written (site 1)(site 1)(site 2)(site 2).
#' Therefore, a genotype of WWHR would be completely wild-type at the first locus,
#' and at the second locus a compound heterozygote with the SEM allele and a
#' costly resistance allele.
#'
#' "V" and "E" alleles are simply minor alleles with the same protein sequence as the
#' major allele at their respective loci, "W". There is the possibility for allelic
#' conversion of the "V" allele into the "W" allele by mechanisms such as MMR.
#'
#' This drive has male and female specific GD and SEM parameters, as well as
#' maternal deposition for both SEM and GD elements. There are no dosage effects
#' modeled (i.e., having two GD alleles increasing or decreasing the GD rates).
#'
#'
#'
#'
#'
#'
#'
#'
#' @param pF Rate of cleavage during GD process in females
#' @param qF Rate of HDR during GD process in females
#' @param rF Rate of in-frame resistance generation during GD process in females
#'
#' @param aF Rate of cleavage during SEM process on GD allele in females
#' @param bF Rate of SSA during SEM process on GD allele in females
#' @param cF Rate of "V" allele formation from SSA during SEM process on GD allele in females
#'
#' @param xF Rate of cleavage during GD process on SEM allele in females
#' @param yF Rate of SSA during GD process on SEM allele in females
#'
#'
#' @param pM Rate of cleavage during GD process in males
#' @param qM Rate of HDR during GD process in males
#' @param rM Rate of in-frame resistance generation during GD process in males
#'
#' @param aM Rate of cleavage during SEM process on GD allele in males
#' @param bM Rate of SSA during SEM process on GD allele in males
#' @param cM Rate of "V" allele formation from SSA during SEM process on GD allele in males
#'
#' @param xM Rate of cleavage during GD process on SEM allele in males
#' @param yM Rate of SSA during GD process on SEM allele in males
#'
#' @param mmrF Rate of MMR in females, driving allelic conversion of "V" into "W"
#' @param mmrM Rate of MMR in males, driving allelic conversion of "V" into "W"
#'
#'
#'
#'
#'
#' @param pDep Rate of cleavage during maternal deposition
#' @param qDep Rate of HDR during maternal deposition
#' @param rDep Rate of in-frame resistance generation during maternal deposition
#'
#'
#'
#'
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
#' @return Named list containing the inheritance cube, transition matrix, genotypes,
#' wild-type allele, and all genotype-specific parameters.
#' @export
cubeSEMtrans <- function(pF=1, qF=1, rF=0,
                    aF=1, bF=1, cF=1,
                    pM=pF, qM=qF, rM=rF,
                    aM=aF, bM=bF, cM=cF,
                    mmrF=0, mmrM=0,
                    pDep=0, qDep=0, rDep=0,
                    eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){




  ## safety checks
  inputVec <- c(pF,qF,rF, aF,bF,cF, pM,qM,rM, aM,bM,cM, mmrF,mmrM, pDep,qDep,rDep)
  if(any(inputVec>1) || any(inputVec<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }








  # # Testing Probs
  # testVec <- runif(n = 19, min = 0, max = 1)
  #
  # mmrF <- testVec[1]; mmrM <- testVec[2]
  # pF <- testVec[3]; qF <- testVec[4]; rF <- testVec[5];
  # aF <- testVec[6]; bF <- testVec[7]; cF <- testVec[8]; dF <- testVec[9];
  # pM <- testVec[10]; qM <- testVec[11]; rM <- testVec[12];
  # aM <- testVec[13]; bM <- testVec[14]; cM <- testVec[15]; dM <- testVec[16];
  # pDep <- testVec[17]; qDep <- testVec[18]; rDep <- testVec[19];











  #############################################################################
  ## generate genotypes
  #############################################################################
  # # List of possible alleles
  # alleles <- list(c('W', 'G', 'U', 'R', 'V', 'S'),
  #                 c('W', 'H', 'R', 'E') )
  # # Generate alleles
  # alleleList <- vector(mode = "list", length = length(alleles))
  # for(i in 1:length(alleles)){
  #   # This generates all combinations of alleles
  #   hold <- as.matrix(x = expand.grid(alleles[[i]], alleles[[i]],
  #                                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  #   # Sort both alleles to remove uniques
  #   #  The alleles aren't unique, and can occur in any order, so we sort them
  #   #  alphabetically in order to remove duplicate genotypes
  #   sortedAlleles <- apply(X = hold, MARGIN = 1, FUN = sort)
  #   # Paste together and get unique ones only
  #   alleleList[[i]] <- unique(do.call(what = paste0, args = list(sortedAlleles[1, ], sortedAlleles[2, ])))
  # }
  # # Get all combinations of both loci
  # openAlleles <- expand.grid(alleleList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # # Paste them together. This is the final list of genotypes
  # genotypes <- file.path(openAlleles[,1], openAlleles[,2], fsep = "")

  genotypes <- c('WWWW','GWWW','UWWW','RWWW','VWWW','SWWW','GGWW','GUWW','GRWW','GVWW',
                 'GSWW','UUWW','RUWW','UVWW','SUWW','RRWW','RVWW','RSWW','VVWW','SVWW',
                 'SSWW','WWHW','GWHW','UWHW','RWHW','VWHW','SWHW','GGHW','GUHW','GRHW',
                 'GVHW','GSHW','UUHW','RUHW','UVHW','SUHW','RRHW','RVHW','RSHW','VVHW',
                 'SVHW','SSHW','WWRW','GWRW','UWRW','RWRW','VWRW','SWRW','GGRW','GURW',
                 'GRRW','GVRW','GSRW','UURW','RURW','UVRW','SURW','RRRW','RVRW','RSRW',
                 'VVRW','SVRW','SSRW','WWEW','GWEW','UWEW','RWEW','VWEW','SWEW','GGEW',
                 'GUEW','GREW','GVEW','GSEW','UUEW','RUEW','UVEW','SUEW','RREW','RVEW',
                 'RSEW','VVEW','SVEW','SSEW','WWHH','GWHH','UWHH','RWHH','VWHH','SWHH',
                 'GGHH','GUHH','GRHH','GVHH','GSHH','UUHH','RUHH','UVHH','SUHH','RRHH',
                 'RVHH','RSHH','VVHH','SVHH','SSHH','WWHR','GWHR','UWHR','RWHR','VWHR',
                 'SWHR','GGHR','GUHR','GRHR','GVHR','GSHR','UUHR','RUHR','UVHR','SUHR',
                 'RRHR','RVHR','RSHR','VVHR','SVHR','SSHR','WWEH','GWEH','UWEH','RWEH',
                 'VWEH','SWEH','GGEH','GUEH','GREH','GVEH','GSEH','UUEH','RUEH','UVEH',
                 'SUEH','RREH','RVEH','RSEH','VVEH','SVEH','SSEH','WWRR','GWRR','UWRR',
                 'RWRR','VWRR','SWRR','GGRR','GURR','GRRR','GVRR','GSRR','UURR','RURR',
                 'UVRR','SURR','RRRR','RVRR','RSRR','VVRR','SVRR','SSRR','WWER','GWER',
                 'UWER','RWER','VWER','SWER','GGER','GUER','GRER','GVER','GSER','UUER',
                 'RUER','UVER','SUER','RRER','RVER','RSER','VVER','SVER','SSER','WWEE',
                 'GWEE','UWEE','RWEE','VWEE','SWEE','GGEE','GUEE','GREE','GVEE','GSEE',
                 'UUEE','RUEE','UVEE','SUEE','RREE','RVEE','RSEE','VVEE','SVEE','SSEE')


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
  #  deposition. Deposition of GD and SEM are indepdent, so they occur "at the same time".

  ##########
  # Female
  ##########
  mendF <- list('one' = list('W' = c('W'=1),
                             'U' = c('U'=1),
                             'R' = c('R'=1),
                             'V' = c('W'=mmrF, 'V'=1-mmrF)),
                'two' = list('W'= c('W'=1),
                             'H'= c('H'=1),
                             'R'= c('R'=1),
                             'E'= c('R'=1)) )

  # assume that G and S are equally competent at homing
  gdF <- list('W' = list('G' = c('W'=1-pF,
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
              'S' = c('S'=1))

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "W" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  semF <- list('one' = list('W' = c('W'=1),
                            'G' = c('G'=1-aF,
                                    'V'=aF*bF*cF,
                                    'W'=aF*bF*(1-cF),
                                    'S'=aF*(1-bF) ),
                            'U' = c('U'=1),
                            'R' = c('R'=1),
                            'V' = c('V'=1),
                            'S' = c('S'=1)),
               'two' = list('W' = c('W'=1),
                            'H' = c('H'=1-xF,
                                    'R'=xF*(1-yF),
                                    'E'=xf*yF),
                            'R' = c('R'=1),
                            'E' = c('E'=1)) )


  ##########
  # Male
  ##########
  mendM <- list('one' = list('W' = c('W'=1),
                             'U' = c('U'=1),
                             'R' = c('R'=1),
                             'V' = c('W'=mmrM, 'V'=1-mmrM)),
                'two' = list('W'= c('W'=1),
                             'H'= c('H'=1),
                             'R'= c('R'=1),
                             'E'= c('R'=1)) )

  # assume that G and S are equally competent at homing
  gdM <- list('W' = list('G' = c('W'=1-pM,
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
              'S' = c('S'=1))

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "W" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  semM <- list('one' = list('W' = c('W'=1),
                            'G' = c('G'=1-aM,
                                    'V'=aM*bM*cM,
                                    'W'=aM*bM*(1-cM),
                                    'S'=aM*(1-bM) ),
                            'U' = c('U'=1),
                            'R' = c('R'=1),
                            'V' = c('V'=1),
                            'S' = c('S'=1)),
               'two' = list('W' = c('W'=1),
                            'H' = c('H'=1-xM,
                                    'R'=xM*(1-yM),
                                    'E'=xM*yM),
                            'R' = c('R'=1),
                            'E' = c('E'=1)) )
















  list(c('W', 'G', 'U', 'R', 'V', 'S'),
                   c('W', 'H', 'R', 'E') )


  # only doing GD deposition - can assume that SEM is controlled and not under
  #  germline expression.
  #  GD deposition only impacts the 'W' allele, so it's the only one we need to
  #  worry about
  # Since we don't know what happened in the female before we got here, we can
  #  have nearly any allele to repair against. The only one we shouldn't see is
  #  "H", so I'll leave it out, and if we do see it, figure out why.
  # Since this is HDR dependent, I give the "V" one an opportunity for increased
  #  MMR, thereby double converting new "V" alleles on into "W" alleles
  dep <- list('W' = c('W'=1-pDep + pDep*qDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'G' = c('W'=1-pDep,
                      'G'=pDep*qDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'U' = c('W'=1-pDep,
                      'U'=pDep*qDep + pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'R' = c('W'=1-pDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*qDep + pDep*(1-qDep)*(1-rDep)),
              'V' = c('W'=1-pDep + pDep*qDep*mmrM,
                      'V'=pDep*qDep*(1-mmrM),
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)),
              'S' = c('W'=1-pDep,
                      'S'=pDep*qDep,
                      'U'=pDep*(1-qDep)*rDep,
                      'R'=pDep*(1-qDep)*(1-rDep)) )
















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
    gdScoreF <- any("G" == fSplit) || any("S" == fSplit)


    ##########
    # Female Alleles
    ##########
    # First step
    if(!gdScoreF){
      # Mendelian both loci
      fAllele1 <- c(mendF[['one']][[ fSplit[1] ]], mendF[['one']][[ fSplit[2] ]])

    } else {
      # GD at first locus
      # Use mendelian at second locus, it will get caught under SEM statement

      # check what the second allele is
      #  If "W", it has a different shape, which depends on which GD allele is
      #  present. Since "W" is only in the second place, we can check it easily.
      if(fSplit[[2]]=='W'){
        # There is a "W" in the second spot
        fAllele1 <- c(gdF[[ fSplit[1] ]], gdF[['W']][[ fSplit[1] ]])
      } else {
        # There is not a "W" in the second spot
        fAllele1 <- c(gdF[[ fSplit[1] ]], gdF[[ fSplit[2] ]])
      }

    } # end first step, Mendelian vs GD

    # handle Mendelian at the second locus
    #  if there is homing or SEM, it will get caught next
    fAllele2 <- c(mendF[['two']][[ fSplit[3] ]], mendF[['two']][[ fSplit[4] ]])

    # Second Step
    if(semScoreF && gdScoreF){
      # SEM
      # Here, we take each allele from above, and we target it with the SEM construct
      #  This is done by multiplying the breakdown of GD results by possible SEM results.
      # Finally, unlist so it is the same shape as the previous two.
      # "Map()" is faster, but I can't figure out the output naming.
      fAllele1 <- unlist(x = mapply(FUN = "*", fAllele1, semF[['one']][names(fAllele1)],
                                   SIMPLIFY = FALSE, USE.NAMES = FALSE),
                        recursive = TRUE, use.names = TRUE)
      fAllele2 <- unlist(x = mapply(FUN = "*", fAllele2, semF[['two']][names(fAllele2)],
                                   SIMPLIFY = FALSE, USE.NAMES = FALSE),
                        recursive = TRUE, use.names = TRUE)
    }

    # Aggregate duplicate alleles
    #  This is mostly useful for SEM, and I could have included it in that step.
    #  However, it will also be often useful in the GD step (but not always)
    #  and is useful in one case in the Mendelian step
    fAlleleReduc1 <- vapply(X = unique(names(fAllele1)),
                           FUN = function(x){sum(fAllele1[names(fAllele1)==x])},
                           FUN.VALUE = numeric(length = 1L))
    fAlleleReduc2 <- vapply(X = unique(names(fAllele2)),
                           FUN = function(x){sum(fAllele2[names(fAllele2)==x])},
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
      gdScoreM <- any("G" == mSplit) || any("S" == mSplit)





list(c('W', 'G', 'U', 'R', 'V', 'S'),
                   c('W', 'H', 'R', 'E') )

      depScoreM <- any("W" == mSplit)







      ##########
      # Male Alleles
      ##########
      # First step
      if(!gdScoreM){
        # Mendelian
        mAllele1 <- c(mendM[['one']][[ mSplit[1] ]], mendM[['one']][[ mSplit[2] ]])

      } else {
        # GD

        # check what the second allele is
        #  If "W", it has a different shape, which depends on which GD allele is
        #  present. Since "W" is only in the second place, we can check it easily.
        if(mSplit[[2]]=='W'){
          # There is a "W" in the second spot
          mAllele1 <- c(gdM[[ mSplit[1] ]], gdM[['W']][[ mSplit[1] ]])
        } else {
          # There is not a "W" in the second spot
          mAllele1 <- c(gdM[[ mSplit[1] ]], gdM[[ mSplit[2] ]])
        }

      } # end first step, Mendelian vs GD

      # handle Mendelian at the second locus
      #  if there is homing or SEM, it will get caught next
      mAllele2 <- c(mendM[['two']][[ mSplit[3] ]], mendM[['two']][[ mSplit[4] ]])

      # Second Step
      if(semScoreM && gdScoreM){
        # SEM
        # Here, we take each allele from above, and we target it with the SEM construct
        #  This is done by multiplying the breakdown of GD results by possible SEM results.
        # Finally, unlist so it is the same shape as the previous two.
        # "Map()" is faster, but I can't figure out the output naming.
        mAllele1 <- unlist(x = mapply(FUN = "*", mAllele1, semM[['one']][names(mAllele1)],
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE),
                          recursive = TRUE, use.names = TRUE)
        mAllele2 <- unlist(x = mapply(FUN = "*", mAllele2, semM[['two']][names(mAllele2)],
                                     SIMPLIFY = FALSE, USE.NAMES = FALSE),
                          recursive = TRUE, use.names = TRUE)
      }

      # Aggregate duplicate alleles
      #  This is mostly useful for SEM, and I could have included it in that step.
      #  However, it will also be often useful in the GD step (but not always)
      #  and is useful in one case in the Mendelian step
      mAlleleReduc1 <- vapply(X = unique(names(mAllele1)),
                             FUN = function(x){sum(mAllele1[names(mAllele1)==x])},
                             FUN.VALUE = numeric(length = 1L))
      mAlleleReduc2 <- vapply(X = unique(names(mAllele2)),
                             FUN = function(x){sum(mAllele2[names(mAllele2)==x])},
                             FUN.VALUE = numeric(length = 1L))








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
          if(mARNames[eleM] == 'W'){
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
    wildType = "WW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "GG"
  ))

}
