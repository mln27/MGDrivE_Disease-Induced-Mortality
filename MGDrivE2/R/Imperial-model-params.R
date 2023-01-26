#------------------------------------------------
#' Model Parameter List Creation
#'
#' \code{model_param_list_create} creates list of model parameters to be used
#' within \code{equilibrium_init_create}
#'
#' @param eta Death rate for expoential population distribtuion, i.e. 1/Mean Population Age. Default = 0.0001305
#' @param rho Age-dependent biting parameter. Default = 0.85
#' @param a0 Age-dependent biting parameter. Default = 2920
#' @param sigma2 Variance of the log heterogeneity in biting rates. Default = 1.67
#' @param max_age Maximum age in days. Default = 100*365
#' @param rA Rate of leaving asymptomatic infection. Default = 0.00512821
#' @param rT Rate of leaving treatment. Default = 0.2
#' @param rD Rate of leaving clinical disease. Default = 0.2
#' @param rU Rate of recovering from subpatent infection. Default = 0.00906627
#' @param rP Rate of leaving prophylaxis. Default = 0.06666667
#' @param dE Latent period of human infection. Default = 12
#' @param delayGam Lag from parasites to infectious gametocytes. Default = 12.5
#' @param cD Untreated disease contribution to infectiousness. Default = 0.0676909
#' @param cT Treated disease contribution to infectiousness. Default =   0.322 * cD
#' @param cU Subpatent disease contribution to infectiousness. Default = 0.006203
#' @param gamma1 Parameter for infectiousness of state A. Default = 1.82425
#' @param d1 Minimum probability due to maximum immunity. Default = 0.160527
#' @param dID Inverse of decay rate. Default = 3650
#' @param ID0 Scale parameter. Default = 1.577533
#' @param kD Shape parameter. Default = 0.476614
#' @param uD Duration in which immunity is not boosted. Default = 9.44512
#' @param aD Scale parameter relating age to immunity. Default = 8001.99
#' @param fD0 Time-scale at which immunity changes with age. Default = 0.007055
#' @param gammaD Shape parameter relating age to immunity. Default = 4.8183
#' @param alphaA PCR detection probability parameters state A. Default = 0.757
#' @param alphaU PCR detection probability parameters state U. Default = 0.186
#' @param b0 Maximum probability due to no immunity. Default = 0.590076
#' @param b1 Maximum relative reduction due to immunity. Default = 0.5
#' @param dB Inverse of decay rate. Default = 3650
#' @param IB0 Scale parameter. Default = 43.8787
#' @param kB Shape parameter. Default = 2.15506
#' @param uB Duration in which immunity is not boosted. Default = 7.19919
#' @param phi0 Maximum probability due to no immunity. Default = 0.791666
#' @param phi1 Maximum relative reduction due to immunity. Default = 0.000737
#' @param dCA Inverse of decay rate. Default = 10950
#' @param IC0 Scale parameter. Default = 18.02366
#' @param kC Shape parameter. Default = 2.36949
#' @param uCA Duration in which immunity is not boosted. Default = 6.06349
#' @param PM New-born immunity relative to mother’s. Default = 0.774368
#' @param dCM Inverse of decay rate of maternal immunity. Default = 67.6952
#' @param tau1 Duration of host seeking, assumed to be constant between species. Default = 0.69
#' @param tau2 Duration of mosquito resting after feed. Default = 2.31
#' @param muF Daily mortality of adult mosquitos. Default = 0.132
#' @param Q0 Anthrophagy probability. Default = 0.92
#' @param nEIP Number of Erlang-distributed EIP compartments. Default = 6
#' @param qEIP Inverse of the mean duration of the EIP. Default = 1/10 (days)
#' @param DY number of days in a year
#' @param thetaB proportion of bites on a person in bed
#' @param thetaI proportion of bites on a person outdoors
#' @param r_llin probability of repeating a feeding attempt due to LLINs
#' @param s_llin probability of feeding and surviving in presence of LLINs
#' @param r_irs probability of repeating a feeding attempt due to IRS
#' @param s_irs probability of feeding and surviving in presence of IRS
#' @param qE mosquito egg lifecycle parameter
#' @param  nE mosquito egg lifecycle parameter
#' @param  qL mosquito larval lifecycle parameter
#' @param  nL mosquito larval lifecycle parameter
#' @param  qP mosquito pupae lifecycle parameter
#' @param  nP mosquito pupae lifecycle parameter
#' @param  muE death rate of egg stage
#' @param  muL death rate of larval stage
#' @param  muP death rate of pupae stage
#' @param  muM death rate of male adult stage
#' @param  eps eggs laid per day
#' @param  nu mosquito lifecycle parameter
#' @param  NH number of humans
#' @param ... Any other parameters needed for non-standard model. If they share the same name
#' as any of the defined parameters \code{model_param_list_create} will stop. You can either write
#' any extra parameters you like individually, e.g. model_param_list_create(extra1 = 1, extra2 = 2)
#' and these parameteres will appear appended to the returned list, or you can pass explicitly
#' the ellipsis argument as a list created before, e.g. model_param_list_create(...=list(extra1 = 1, extra2 = 2))
#'
#'
#'
#' This function creates all of the necessary parameters for the Imperial model. Parameters furnished by MGDrivE will be
#' removed from this function. Adapted from: https://github.com/mrc-ide/deterministic-malaria-model/blob/master/R/model_parameters.R
#'
#' @export
imperial_model_param_list_create <- function(
  # age, heterogeneity in exposure,
  eta = 1/(21*365),
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.67,
  max_age = 100*365,
  #  rate of leaving infection states
  rA = 1/195,
  rT = 0.2,
  rD = 0.2,
  rU = 1/110.299,
  rP = 1/15,
  #  human latent period and time lag from asexual parasites to
  dE  = 12,
  delayGam = 12.5,
  # human infectiousness to mosquitoes
  cD  = 0.0676909,
  cT  =  0.322 * cD,
  cU  = 0.006203,
  gamma1  = 1.82425,
  #  Immunity reducing probability of detection
  d1 = 0.160527,
  dID = 3650,
  ID0 = 1.577533,
  kD = 0.476614,
  uD = 9.44512,
  aD = 8001.99,
  fD0 = 0.007055,
  gammaD = 4.8183,
  alphaA = 0.75735,
  alphaU = 0.185624,
  # Immunity reducing probability of infection
  b0 = 0.590076,
  b1 = 0.5,
  dB = 3650,
  IB0 = 43.8787,
  kB = 2.15506,
  uB = 7.19919,
  # Immunity reducing probability of clinical disease
  phi0 = 0.791666,
  phi1 = 0.000737,
  dCA = 10950,
  IC0 = 18.02366,
  kC = 2.36949,
  uCA = 6.06349,
  PM = 0.774368,
  dCM = 67.6952,
  # entomological parameters
  tau1 = 0.69,
  tau2 = 2.31,
  muF = 0.132, # female mortality
  nEIP = 3, # number of EIP compartments
  qEIP = 1/10, # inverse of mean duration of EIP
  Q0 = 0.92,
  DY = 365,
  thetaB = 0.89, # proportion of bites on a person in bed
  thetaI = 0.97, # proportion of bites on a person outdoors
  r_llin = 0.56, # probability of repeating a feeding attempt due to LLINs
  s_llin = 0.03, # probability of feeding and surviving in presence of LLINs
  r_irs = 0.60, # probability of repeating a feeding attempt due to IRS
  s_irs = 0, # probability of feeding and surviving in presence of IRS
  # lifecycle parameters
  qE = 1/4,
  nE = 2,
  qL = 1/3,
  nL = 3,
  qP = 1/6,
  nP = 2,
  muE = 0.05,
  muL = 0.15,
  muP = 0.05,
  muM = 0.132,
  eps = 58.9,
  nu = 1/(4/24),
  # epidemiological parameters
  NH = 1000,
  ...

){
  # set up param list
  mp_list <- list()

  # catch extra params and place in list
  extra_param_list <- list(...)
  if(length(extra_param_list)>0){
    if(is.list(extra_param_list[[1]])){
  extra_param_list <- extra_param_list[[1]]
    }
  }

  ## DEFAULT PARAMS

  # duration of year
  mp_list$DY <- DY

  # age, heterogeneity in exposure
  mp_list$eta <- eta
  mp_list$rho <- rho
  mp_list$a0 <- a0
  mp_list$sigma2 <- sigma2
  mp_list$max_age <- max_age

  # rate of leaving infection states
  mp_list$rA <- rA
  mp_list$rT <- rT
  mp_list$rD <- rD
  mp_list$rU <- rU
  mp_list$rP <- rP

  # human latent period and time lag from asexual parasites to
  # infectiousness
  mp_list$dE <- dE
  mp_list$delayGam <- delayGam

  # infectiousness to mosquitoes
  mp_list$cD <- cD
  mp_list$cT <- cT
  mp_list$cU <- cU
  mp_list$gamma1 <- gamma1

  # Immunity reducing probability of detection
  mp_list$d1 <- d1
  mp_list$dID <- dID
  mp_list$ID0 <- ID0
  mp_list$kD <- kD
  mp_list$uD <- uD
  mp_list$aD <- aD
  mp_list$fD0 <- fD0
  mp_list$gammaD <- gammaD

  # PCR prevalence parameters
  mp_list$alphaA <- alphaA
  mp_list$alphaU <- alphaU

  # anti-infection immunity
  mp_list$b0 <- b0
  mp_list$b1 <- b1
  mp_list$dB <- dB
  mp_list$IB0 <- IB0
  mp_list$kB <- kB
  mp_list$uB <- uB

  # clinical immunity
  mp_list$phi0 <- phi0
  mp_list$phi1 <- phi1
  mp_list$dCA <- dCA
  mp_list$IC0 <- IC0
  mp_list$kC <- kC
  mp_list$uCA <- uCA
  mp_list$PM <- PM
  mp_list$dCM <- dCM

  # entomological parameters

    # lifecycle parameters
  mp_list$qE <- qE
  mp_list$nE <- nE
  mp_list$qL <- qL
  mp_list$nL <- nL
  mp_list$qP <- qP
  mp_list$nP <- nP
  mp_list$muE <- muE
  mp_list$muL <- muL
  mp_list$muP <- muP
  mp_list$muF <- muF
  mp_list$muM <- muM
  mp_list$eps <- eps
  mp_list$nu <- nu

  # epidemiological parameters
  mp_list$NH <- NH
  mp_list$nEIP <- nEIP
  mp_list$qEIP <- qEIP
  mp_list$tau1 <- tau1
  mp_list$tau2 <- tau2
  mp_list$muF <- muF
  mp_list$Q0 <- Q0
  mp_list$fv0 <- 1 / (tau1 + tau2)
  mp_list$av0 <- Q0 * mp_list$fv0 # daily feeeding rate on humans
  # in Erlang distributed EIP, survival probability is proportional to number of
  # compartments, death rate, and mean duration in each compartment
  mp_list$Surv0 <- (mp_list$nEIP * mp_list$qEIP/ (mp_list$nEIP * mp_list$qEIP + mp_list$muF))^mp_list$nEIP
  mp_list$beta <- (eps*muF)/(exp(muF/mp_list$fv0)-1)


  # intervention related parameters
  mp_list$thetaB <- thetaB
  mp_list$thetaI <- thetaI
  mp_list$r_llin <- r_llin
  mp_list$s_llin <- s_llin
  mp_list$r_irs <- r_irs
  mp_list$s_irs <- s_irs

  return(append(mp_list,extra_param_list))
}
