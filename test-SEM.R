###############################################################################
#   _    _____       _____ ________  ___
#  | |  / /__ \     / ___// ____/  |/  /
#  | | / /__/ /_____\__ \/ __/ / /|_/ /
#  | |/ // __/_____/__/ / /___/ /  / /
#  |___//____/    /____/_____/_/  /_/
# .............................................................................
# V2 SEM Test Script
# Ndeffo Mbah Lab
# Jared Bennett, jared_bennett@berkeley.edu
#
# 20220328
# Initial setup.
#
#
###############################################################################
###############################################################################
####################
## Setup Packages
####################
rm(list=ls());gc()
set.seed(10)
library(MGDrivE2)
library(ggplot2)


####################
## Parameters
####################
# number of adult female mosquitoes
NF <- 500
tmax <- 365
dt <- 1
dt_stoch <- 0.1
nReps <- 3

# entomological parameters
theta <- list(
  qE = 1/4,
  nE = 2,
  qL = 1/3,
  nL = 3,
  qP = 1/6,
  nP = 2,
  muE = 0.05,
  muL = 0.15,
  muP = 0.05,
  muF = 0.09,
  muM = 0.09,
  beta = 16,
  nu = 1/(4/24)
)

####################
## Inheritance
####################
# no MMR of V to W
# no maternal deposition
# equal parameters between males and females
cube <- MGDrivE2::cubeSEM(pF = 0.8, qF = 0.4, rF = 0.3,
                          aF = 0.95, bF = 0.95, cF = 0.6, dF = 0.001)




####################
## Releases
####################
# gene drive releases
gdEvents <- data.frame("var" = paste0("F_", cube$releaseType, "_", cube$wildType),
                       "time" = seq(from = 20, length.out = 5, by = 10),
                       "value" = 50,
                       "method" = "add",
                       stringsAsFactors = FALSE)

# SEM "releases"





cube$genotypesID








####################
## Setup Petri Net
####################
# Places
SPN_P <- spn_P_lifecycle_node(params = theta,cube = cube)

# Transitions
SPN_T <- spn_T_lifecycle_node(spn_P = SPN_P,params = theta,cube = cube)

# Stoichiometry matrix
S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)


####################
## Equilibrium and Hazards
####################
# calculate equilibrium and setup initial conditions
initialCons <- equilibrium_lifeycle(params = theta, NF = NF, phi = 0.5,
                                    log = TRUE, spn_P = SPN_P, cube = cube)

# approximate hazards for continous approximation
approx_hazards <- spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube,
                              params = initialCons$params, type = "life",
                              log = TRUE, exact = FALSE, tol = 1e-8,
                              verbose = FALSE)


SPN_T$T[[5569]]

length(SPN_T$T)
length(SPN_T$v)

hold <- unlist(lapply(SPN_T$T, '[[', "class"))

which(is.null(hold))





# exact hazards for integer-valued state space
exact_hazards <- spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube,
                             params = initialCons$params, type = "life",
                             log = TRUE, exact = TRUE, tol = NaN,
                             verbose = FALSE)


####################
## Simulate
####################
# deterministic
ODE_out <- sim_trajectory_R(x0 = initialCons$M0, tmax = tmax, dt = dt, S = S,
                            hazards = approx_hazards, sampler = "ode", method = "lsoda",
                            events = events, verbose = FALSE)

# tau sampling
PTS_out <- sim_trajectory_R(x0 = initialCons$M0, tmax = tmax, dt = dt, num_reps = nReps,
                            dt_stoch = dt_stoch, S = S, hazards = exact_hazards,
                            sampler = "tau", events = events, verbose = FALSE)


####################
## Analyze and Plot
####################
##########
# Deterministic
##########
# summarize females by genotype
ODE_out_f <- summarize_females(out = ODE_out$state, spn_P = SPN_P)

# summarize males by genotype
ODE_out_m <- summarize_males(out = ODE_out$state)

# add sex for plotting
ODE_out_f$sex <- "Female"
ODE_out_m$sex <- "Male"

# plot
ggplot(data = rbind(ODE_out_f, ODE_out_m)) +
  geom_line(aes(x = time, y = value, color = genotype)) +
  facet_wrap(facets = vars(sex), scales = "fixed") +
  theme_bw() +
  ggtitle("SPN: ODE Solution")


##########
# Stochastic
##########




# summarize females/males
PTS_out_f <- summarize_females(out = PTS_out$state, spn_P = SPN_P)
PTS_out_m <- summarize_males(out = PTS_out$state)

# add sex for plotting
PTS_out_f$sex <- "Female"
PTS_out_m$sex <- "Male"

# plot adults
ggplot(data = rbind(PTS_out_f, PTS_out_m)) +
  geom_line(aes(x = time, y = value, color = genotype)) +
  facet_wrap(facets = vars(sex), scales = "fixed") +
  theme_bw() +
  ggtitle("SPN: Tau-leaping Approximation")



