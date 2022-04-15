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
#   April 2022
#
###############################################################################

#' Inheritance Cube: X-Linked ReMEDE (REpeat Mediated Excision of a Drive Element) in Trans
#'
#' This is the second ReMEDE system put forth by \href{GETLINKWHENEXISTS!!!}{Chennuri and Myles}.
#' It's another realization of \href{https://doi.org/10.1098/rstb.2019.0804}{biodegradable gene drives},
#' where the pieces exist as independently-recombining constructs that only act
#' when in the presence of each other.
#'
#' This construct consists of one X-linked target site and one autosomal target site.
#' The first site, on the X chromosome, carries
#' gRNA and Cas9 for homing and removal of the SEM site. The second site, on
#' an autosome, consists
#' of an endonuclease targeting the Cas9 site. Both sites contain direct repeats
#' to induce SSA upon chromosome damage. There are 7 possible alleles at the first site:
#' \itemize{
#'  \item X: Wild-type X allele
#'  \item G: Active GD (gene drive), with targetable SEM (self-elimination mechanism) site
#'  \item U: Low-cost resistant allele, non-targetable by GD or SEM (R1 in other literature)
#'  \item R: High-cost resistant allele, non-targetable by GD or SEM (R2 in other literature)
#'  \item V: "wild-type", product of SEM, non-targetable by GD or SEM
#'  \item S: Active GD with non-targetable SEM site
#'  \item Y: Y-chromosome, non-targetable
#' }
#'
#' There are 4 possible alleles at the autosomal second site:
#' \itemize{
#'  \item W: Wild-type allele
#'  \item H: Active SEM element
#'  \item R: High-cost resistant allele, non-targetable by GD (R2 in other literature)
#'  \item E: "wild-type", produce of GD, non-targetbale by GD
#' }
#'
#' This provides a total genotype count of 270 in the cube, however, only 210 viable
#' female genotypes and 60 viable male genotypes. Genotypes are written (site 1)(site 1)(site 2)(site 2),
#' where "site 1" is X-linked and "site 2" is autosomally linked.
#' Therefore, a genotype of XXHR would be female and completely wild-type at the first locus,
#' and at the second locus a compound heterozygote with the SEM allele and a
#' costly resistance allele.
#'
#' "V" and "E" alleles are simply minor alleles with the same protein sequence as the
#' major allele at their respective loci, "W" or "X". There is the possibility for allelic
#' conversion of the "V" allele into the "X" allele by mechanisms such as MMR.
#'
#' This drive has male and female specific GD and SEM parameters, as well as
#' maternal deposition for both SEM and GD elements. There are no dosage effects
#' modeled (i.e., having two GD alleles increasing or decreasing the GD rates).
#'
#' Gene-drive and SEM parameters (p*, q*, r*, a*, b*, c*, x*, y*, mmr*) are
#' all rates and values must fall in the range of [0, 1], inclusive.
#'
#' Population parameters (phi, xiF, xiM) are either NULL, for default values, or
#' named vectors, where the names must match the names in the inheritance pattern.
#' These parameters are also rates and fall in the range of [0, 1], inclusive.
#'
#' Omega and s may also be NULL or named vectors, but are weights, and thus fall
#' in the range [0, inf).
#'
#' Eta, the male mating weights, may be NULL for default values, or applied as a
#' list of vectors. If the vectors are length 2 (e.g., c(maleGeno, matingWeight)),
#' then that weight is applied to all female genotypes for that specific male genotype.
#' Alternatively, vectors of length 3 may be supplied (e.g., c(femaleGeno, maleGeno, matingWeight)).
#' In this case, the mating weight only applies to the specific female X male mating
#' event. The two lists may NOT be mixed - the list must hold all length 2 or all
#' length 3 vectors.
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
#' @param pDep Rate of cleavage during maternal deposition into W allele from GD
#' @param qDep Rate of HDR during maternal deposition into W allele from GD
#' @param rDep Rate of in-frame resistance generation during maternal deposition into W allele from GD
#'
#' @param aDep Rate of cleavage during maternal deposition into G allele from SEM
#' @param bDep Rate of SSA during maternal deposition into G allele from SEM
#' @param cDep Rate of MMR, converting G into V, from maternal deposition from SEM
#'
#' @param xDep Rate of cleavage of SEM allele from GD during maternal deposition
#' @param yDep Rate of SSA of SEM allele from GD during maternal deposition
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
cubeSEMtransX_GD <- function(pF=1, qF=1, rF=0,
                             aF=1, bF=1, cF=1,
                             xF=1, yF=1,
                             aM=aF, bM=bF, cM=cF,
                             xM=xF, yM=yF,
                             mmrF=0, mmrM=mmrF,
                             pDep=0, qDep=0, rDep=0,
                             aDep=0, bDep=0, cDep=0,
                             xDep=0, yDep=0,
                             eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  ## safety checks
  inputVec <- c(pF,qF,rF, aF,bF,cF, xF,yF,
                aM,bM,cM, xM,yM,
                mmrF,mmrM,
                pDep,qDep,rDep, aDep,bDep,cDep, xDep,yDep)
  if(any(inputVec>1) || any(inputVec<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }


  # # Testing Probs
  # testVec <- runif(n = 26, min = 0, max = 1)
  #
  # mmrF <- testVec[1]; mmrM <- testVec[2]
  # pF <- testVec[3]; qF <- testVec[4]; rF <- testVec[5];
  # aF <- testVec[6]; bF <- testVec[7]; cF <- testVec[8];
  # xF <- testVec[9]; yF <- testVec[10];
  # #pM <- testVec[11]; qM <- testVec[12]; rM <- testVec[13];
  # aM <- testVec[14]; bM <- testVec[15]; cM <- testVec[16];
  # xM <- testVec[17]; yM <- testVec[18];
  # pDep <- testVec[19]; qDep <- testVec[20]; rDep <- testVec[21];
  # aDep <- testVec[22]; bDep <- testVec[23]; cDep <- testVec[24];
  # xDep <- testVec[25]; yDep <- testVec[26];


  #############################################################################
  ## generate genotypes
  #############################################################################
  # # List of possible alleles
  # alleles <- list(c('X', 'G', 'U', 'R', 'V', 'S', 'Y'),
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
  # # get female genotypes
  # genoF <- grep(pattern = "Y", x = genotypes, value = TRUE, invert = TRUE)
  # # get viable male genotypes
  # genoM <- grep(pattern = "YY",
  #               x = grep(pattern = "Y", x = genotypes, value = TRUE, invert = FALSE),
  #               value = TRUE, invert = TRUE)

  genoF <- c('XXWW','GXWW','UXWW','RXWW','VXWW','SXWW','GGWW','GUWW','GRWW','GVWW','GSWW',
             'UUWW','RUWW','UVWW','SUWW','RRWW','RVWW','RSWW','VVWW','SVWW','SSWW','XXHW',
             'GXHW','UXHW','RXHW','VXHW','SXHW','GGHW','GUHW','GRHW','GVHW','GSHW','UUHW',
             'RUHW','UVHW','SUHW','RRHW','RVHW','RSHW','VVHW','SVHW','SSHW','XXRW','GXRW',
             'UXRW','RXRW','VXRW','SXRW','GGRW','GURW','GRRW','GVRW','GSRW','UURW','RURW',
             'UVRW','SURW','RRRW','RVRW','RSRW','VVRW','SVRW','SSRW','XXEW','GXEW','UXEW',
             'RXEW','VXEW','SXEW','GGEW','GUEW','GREW','GVEW','GSEW','UUEW','RUEW','UVEW',
             'SUEW','RREW','RVEW','RSEW','VVEW','SVEW','SSEW','XXHH','GXHH','UXHH','RXHH',
             'VXHH','SXHH','GGHH','GUHH','GRHH','GVHH','GSHH','UUHH','RUHH','UVHH','SUHH',
             'RRHH','RVHH','RSHH','VVHH','SVHH','SSHH','XXHR','GXHR','UXHR','RXHR','VXHR',
             'SXHR','GGHR','GUHR','GRHR','GVHR','GSHR','UUHR','RUHR','UVHR','SUHR','RRHR',
             'RVHR','RSHR','VVHR','SVHR','SSHR','XXEH','GXEH','UXEH','RXEH','VXEH','SXEH',
             'GGEH','GUEH','GREH','GVEH','GSEH','UUEH','RUEH','UVEH','SUEH','RREH','RVEH',
             'RSEH','VVEH','SVEH','SSEH','XXRR','GXRR','UXRR','RXRR','VXRR','SXRR','GGRR',
             'GURR','GRRR','GVRR','GSRR','UURR','RURR','UVRR','SURR','RRRR','RVRR','RSRR',
             'VVRR','SVRR','SSRR','XXER','GXER','UXER','RXER','VXER','SXER','GGER','GUER',
             'GRER','GVER','GSER','UUER','RUER','UVER','SUER','RRER','RVER','RSER','VVER',
             'SVER','SSER','XXEE','GXEE','UXEE','RXEE','VXEE','SXEE','GGEE','GUEE','GREE',
             'GVEE','GSEE','UUEE','RUEE','UVEE','SUEE','RREE','RVEE','RSEE','VVEE','SVEE',
             'SSEE')
  genoM <- c('XYWW','GYWW','UYWW','RYWW','VYWW','SYWW','XYHW','GYHW','UYHW','RYHW','VYHW',
             'SYHW','XYRW','GYRW','UYRW','RYRW','VYRW','SYRW','XYEW','GYEW','UYEW','RYEW',
             'VYEW','SYEW','XYHH','GYHH','UYHH','RYHH','VYHH','SYHH','XYHR','GYHR','UYHR',
             'RYHR','VYHR','SYHR','XYEH','GYEH','UYEH','RYEH','VYEH','SYEH','XYRR','GYRR',
             'UYRR','RYRR','VYRR','SYRR','XYER','GYER','UYER','RYER','VYER','SYER','XYEE',
             'GYEE','UYEE','RYEE','VYEE','SYEE')
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
  #  deposition. Deposition of GD and SEM are indepdent, so they occur "at the same time".

  ##########
  # Female
  ##########
  mendF <- list('one' = list('X' = c('X'=1),
                             'U' = c('U'=1),
                             'R' = c('R'=1),
                             'V' = c('X'=mmrF, 'V'=1-mmrF)),
                'two' = list('W'= c('W'=1),
                             'H'= c('H'=1),
                             'R'= c('R'=1),
                             'E'= c('R'=1)) )

  # assume that G and S are equally competent at homing
  gdF <- list('X' = list('G' = c('X'=1-pF,
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
              'S' = c('S'=1))

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "W" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  semF <- list('one' = list('X' = c('X'=1),
                            'G' = c('G'=1-aF,
                                    'V'=aF*bF*cF,
                                    'X'=aF*bF*(1-cF),
                                    'S'=aF*(1-bF) ),
                            'U' = c('U'=1),
                            'R' = c('R'=1),
                            'V' = c('V'=1),
                            'S' = c('S'=1)),
               'two' = list('W' = c('W'=1),
                            'H' = c('H'=1-xF,
                                    'R'=xF*(1-yF),
                                    'E'=xF*yF),
                            'R' = c('R'=1),
                            'E' = c('E'=1)) )


  ##########
  # Male
  ##########
  mendM <- list('one' = list('X' = c('X'=1),
                             'G' = c('G'=1),
                             'U' = c('U'=1),
                             'R' = c('R'=1),
                             'V' = c('X'=mmrM, 'V'=1-mmrM),
                             'S' = c('S'=1),
                             'Y' = c('Y'=1)),
                'two' = list('W'= c('W'=1),
                             'H'= c('H'=1),
                             'R'= c('R'=1),
                             'E'= c('R'=1)) )

  # This set assumes gdM has occurred, otherwise some probabilities will be wrong.
  #  "X" isn't targeted a second time
  #  "V" doesn't undergo a second round of MMR
  semM <- list('one' = list('X' = c('X'=1),
                            'G' = c('G'=1-aM,
                                    'V'=aM*bM*cM,
                                    'X'=aM*bM*(1-cM),
                                    'S'=aM*(1-bM) ),
                            'U' = c('U'=1),
                            'R' = c('R'=1),
                            'V' = c('V'=1),
                            'S' = c('S'=1),
                            'Y' = c('Y'=1)),
               'two' = list('W' = c('W'=1),
                            'H' = c('H'=1-xM,
                                    'R'=xM*(1-yM),
                                    'E'=xM*yM),
                            'R' = c('R'=1),
                            'E' = c('E'=1)) )

  # Repair at the first locus depends on the allele and type of deposition
  #  1 - if the mother deposits GD and the male has a 'W' allele there, we have
  #      HDR dependent repair of the 'W' allele.
  #  2 - if the mother deposits SEM and the male has a 'G' allele, we have
  #      SSA dependent repair of the 'G' allele.
  # Thus, this is a 2-level list, the first level dependent on the male allele
  #  ('W' or 'G' currently), and the second level dependent on the maternal allele
  #  (could be any of them).
  # SSA dependent repair does not depend on the maternal allele, but to keep
  #  the same list shape, I have put all of the alleles in. Notice, the outcomes
  #  are all the same in this instance.
  # HDR does depend on the maternal allele for repair, and I have providedthe "V"
  #  allele one opportunity for alleleic conversion directly into the "W" allele.
  dep1 <- list('X' = list('X' = c('X'=1-pDep + pDep*qDep,
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
                                  'R'=pDep*(1-qDep)*(1-rDep)) ),

               'G' = list('X' = c('X'=aDep*bDep*(1-cDep),
                                  'G'=1-aDep,
                                  'V'=aDep*bDep*cDep,
                                  'S'=aDep*(1-bDep)),
                          'G' = c('X'=aDep*bDep*(1-cDep),
                                  'G'=1-aDep,
                                  'V'=aDep*bDep*cDep,
                                  'S'=aDep*(1-bDep)),
                          'U' = c('X'=aDep*bDep*(1-cDep),
                                  'G'=1-aDep,
                                  'V'=aDep*bDep*cDep,
                                  'S'=aDep*(1-bDep)),
                          'R' = c('X'=aDep*bDep*(1-cDep),
                                  'G'=1-aDep,
                                  'V'=aDep*bDep*cDep,
                                  'S'=aDep*(1-bDep)),
                          'V' = c('X'=aDep*bDep*(1-cDep),
                                  'G'=1-aDep,
                                  'V'=aDep*bDep*cDep,
                                  'S'=aDep*(1-bDep)),
                          'S' = c('X'=aDep*bDep*(1-cDep),
                                  'G'=1-aDep,
                                  'V'=aDep*bDep*cDep,
                                  'S'=aDep*(1-bDep)) )
  )

  # repair at the second locus is always SSA, so it actually does NOT depend on the
  #  female allele. However, to keep the structure the same as above, I will
  #  build and use a similarly shaped list.
  # This is a 1-level list, because we only worry about the "H" allele from the
  #  male, and subset the "correct" female allele for pairing.
  dep2 <- list('W' = c('H'=1-xDep,
                       'E'=xDep*yDep,
                       'R'=xDep*(1-yDep)),
               'H' = c('H'=1-xDep,
                       'E'=xDep*yDep,
                       'R'=xDep*(1-yDep)),
               'R' = c('H'=1-xDep,
                       'E'=xDep*yDep,
                       'R'=xDep*(1-yDep)),
               'E' = c('H'=1-xDep,
                       'E'=xDep*yDep,
                       'R'=xDep*(1-yDep))
               )


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
      if(fSplit[[2]]=='X'){
        # There is a "W" in the second spot
        fAllele1 <- c(gdF[[ fSplit[1] ]], gdF[['X']][[ fSplit[1] ]])
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
    for(mi in (1:numGenM + numGenF)){
      # Male loop
      # Split genotype into alleles
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]

      # Score them
      semScoreM <- any("H" == mSplit)
      gScoreM <- any("G" == mSplit) # for deposition later
      gdScoreM <- gScoreM || any("S" == mSplit)
      xScoreM <- any("X" == mSplit) # for deposition later, and X only occurs in first locus
      #  could we just check the first allele? is that the only place where X occurs?


      ##########
      # Male Alleles
      ##########
      # First step
      # Mendelian at the first locus
      mAllele1 <- c(mendM[['one']][[ mSplit[1] ]], mendM[['one']][[ mSplit[2] ]])

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
      #  probability and allele sorting/pasting code into several places:
      #   - no deposition at all
      #   - locus 1, deposition with or without impact
      #   - locus 2, deposition with or without impact
      #  Additionally, we grow a vector in the deposition section, which isn't great
      #  but I don't have a better idea (and it's small, so hopefully not terrible).

      ####################
      # Maternal deposition check
      ####################
      ##########
      # Locus 1
      ##########
      # At locus 1 there can be deposition if:
      #  The mother provides a GD allele, and the father a "W" at the first locus
      #  The mother has an SEM allele, and the father a "G" allele
      # This first check is effectively "are there any "W" or "G" alleles at all"
      #  This means that deposition is possible, not which kind or that it does happen.
      #
      # This may be overbuilt? Am I adding too much for some reason?
      if( (gdScoreF && xScoreM) || (semScoreF && gScoreM) ){
        # there is deposition of some type at locus 1

        # return object
        lProbs1 <- numeric(length = 0L)
        lNames1 <- character(length = 0L)

        # things I reference a lot
        mARNames <- names(mAlleleReduc1)
        fARNames <- names(fAlleleReduc1)

        # loop over male alleles
        for(eleM in 1:length(mAlleleReduc1)){

          # check allele/maternal deposition combination for repair structure
          #  there are 2 possible scenarios, plus the basic "nothing happening" case
          xCheck <- mARNames[eleM] == 'X'
          gCheck <- mARNames[eleM] == 'G'
          if( (gdScoreF && xCheck) || (semScoreF && gCheck) ){
            # grab proper reference for deposition
            #  This has to happen inside the male allele loop because now deposition
            #  depends on precisely which allele we're working on.
            #  We don't want GD deposition from females to pair with a "G" allele
            #  in males, because that shouldn't cause anything to happen.

            # Repair during deposition depends on the maternal allele inherited,
            #  so now we loop over all female alleles, and use that to define
            #  the repair process.
            # The we process everything and append to geno/probs vectors
            for(eleF in 1:length(fAlleleReduc1)){
              # Get repair distribution
              repVec <- dep1[[ mARNames[eleM] ]][[ fARNames[eleF] ]]

              ##########
              # Genotype probs
              ##########
              #  These come from the male allele probs, the female allele probs,
              #  and the distribution of repair products
              lProbs1 <- c(lProbs1, mAlleleReduc1[eleM] * fAlleleReduc1[eleF] * repVec)

              ##########
              # Genotypes
              ##########
              #  Get all combinations of alleles
              oneAGenos <- expand.grid(names(repVec), fARNames[eleF],
                                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
              #  Sort alleles
              oneAGenos <- apply(X = oneAGenos, MARGIN = 1, FUN = sort.int)
              #  paste and append genotypes
              lNames1 <- c(lNames1, file.path(oneAGenos[1,], oneAGenos[2,], fsep = ""))

            } # end female allele loop

          } else {
            # Not a "W" or "G" allele
            # These alleles are NOT impacted by deposition, so their outcome
            #  is simply the probability of the male allele, multiplied by the
            #  probability of each female allele

            ##########
            # Genotype probs
            ##########
            lProbs1 <- c(lProbs1, mAlleleReduc1[eleM] * fAlleleReduc1)

            ##########
            # Genotypes
            ##########
            #  Get all combinations of alleles
            oneAGenos <- expand.grid(mARNames[eleM], fARNames, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
            #  Sort alleles
            oneAGenos <- apply(X = oneAGenos, MARGIN = 1, FUN = sort.int)
            #  paste and append genotypes
            lNames1 <- c(lNames1, file.path(oneAGenos[1,], oneAGenos[2,], fsep = ""))

          } # end male allele check

        } # end male allele loop

      } else {
        # there is not deposition at locus 1
        # Combinations of alleles at locus 1
        holdProbs1 <- expand.grid(fAlleleReduc1, mAlleleReduc1,
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele1 <- expand.grid(names(fAlleleReduc1), names(mAlleleReduc1),
                                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # Sort and collapse
        holdAllele1 <- apply(X = holdAllele1, MARGIN = 1, FUN = sort.int)

        lProbs1 <- holdProbs1[ ,1]*holdProbs1[ ,2]
        lNames1 <- file.path(holdAllele1[1,], holdAllele1[2,], fsep = "")

      } # end locus 1 deposition check


      ##########
      # Locus 2
      ##########
      # At locus 2, the SEM locus, there can be deposition of the mother provides
      #  a GD allele and the father has an SEM allele
      if(gdScoreF && semScoreM){
        # there is maternal deposition from GD at the male SEM locus

        # return objects
        lProbs2 <- numeric(length = 0L)
        lNames2 <- character(length = 0L)

        # things I use a lot
        mARNames <- names(mAlleleReduc2)
        fARNames <- names(fAlleleReduc2)

        # loop over male alleles
        for(eleM in 1:length(mAlleleReduc2)){
          # check for allele type
          #  only "H" alleles are impacted by deposition, so only they have
          #  complicated repair structures
          if(mARNames[eleM] == 'H'){
            # "H" allele

            # Repair during deposition depends on the maternal allele inherited,
            #  so now we loop over all female alleles, and use that to define
            #  the repair process.
            # The we process everything and append to geno/probs vectors
            for(eleF in 1:length(fAlleleReduc2)){
              # Get repair distribution
              repVec <- dep2[[ fARNames[eleF] ]]

              ##########
              # Genotype probs
              ##########
              #  These come from the male allele probs, the female allele probs,
              #  and the distribution of repair products
              lProbs2 <- c(lProbs2, mAlleleReduc2[eleM] * fAlleleReduc2[eleF] * repVec)

              ##########
              # Genotypes
              ##########
              #  Get all combinations of alleles
              oneAGenos <- expand.grid(names(repVec), fARNames[eleF],
                                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
              #  Sort alleles
              oneAGenos <- apply(X = oneAGenos, MARGIN = 1, FUN = sort.int)
              #  paste and append genotypes
              lNames2 <- c(lNames2, file.path(oneAGenos[1,], oneAGenos[2,], fsep = ""))

            } # end female allele loop

          } else {
            # Not an "H" allele
            # These alleles are NOT impacted by deposition, so their outcome
            #  is simply the probability of the male allele, multiplied by the
            #  probability of each female allele

            ##########
            # Genotype probs
            ##########
            lProbs2 <- c(lProbs2, mAlleleReduc2[eleM] * fAlleleReduc2)

            ##########
            # Genotypes
            ##########
            #  Get all combinations of alleles
            oneAGenos <- expand.grid(mARNames[eleM], fARNames, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
            #  Sort alleles
            oneAGenos <- apply(X = oneAGenos, MARGIN = 1, FUN = sort.int)
            #  paste and append genotypes
            lNames2 <- c(lNames2, file.path(oneAGenos[1,], oneAGenos[2,], fsep = ""))

          } # end male allele check

        } # end male allele loop

      } else {
        # At this point, the father must not have an SEM allele, so there's
        #  no deposition at the second locus
        # Combinations of alleles at locus 2
        holdProbs2 <- expand.grid(fAlleleReduc2, mAlleleReduc2,
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele2 <- expand.grid(names(fAlleleReduc2), names(mAlleleReduc2),
                                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # Sort and collapse
        holdAllele2 <- apply(X = holdAllele2, MARGIN = 1, FUN = sort.int)

        lProbs2 <- holdProbs2[ ,1]*holdProbs2[ ,2]
        lNames2 <- file.path(holdAllele2[1,], holdAllele2[2,], fsep = "")

      } # end locus 2 deposition check


      ####################
      # Finish Genotype Distribution
      ####################
      # Expand genotype combinations of locus 1 and locus 2
      holdGNames <- expand.grid(lNames1, lNames2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdGProbs <- expand.grid(lProbs1, lProbs2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # Collapse genotype probs and names
      genoProbs <- holdGProbs[ ,1] * holdGProbs[ ,2]
      genoNames <- file.path(holdGNames[ ,1], holdGNames[ ,2], fsep = "")

      # Aggregate duplicate genotypes
      genoAgg <- vapply(X = unique(genoNames),
                        FUN = function(x){sum(genoProbs[genoNames==x])},
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
    wildType = c('XXWW','XYWW'),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = 'GYWW'
  ))

}
