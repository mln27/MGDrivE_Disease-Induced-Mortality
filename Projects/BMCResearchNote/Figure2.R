###############################################################################
#      ____  __  _________            ___
#     / __ )/  |/  / ____/           |__ \
#    / __  / /|_/ / /      ______    __/ /
#   / /_/ / /  / / /___   /_____/   / __/
#  /_____/_/  /_/\____/            /____/
# .............................................................................
# BMC Research Note - Figure 2
# Ndeffo Mbah Lab
# Jared Bennett, jared_bennett@berkeley.edu
#
# 20220525
#  This script is designed to generate figure 2 of the BCM research note.
#  It is significantly taken from the vignette demonstrating D.I.M. None of the
#  parameters are novel, they all come from V2 vignettes. The point is simply
#  to show that the idea and implementation work, not to demonstrate any kind of
#  scientific point.
#
# 20220526
#  Found a bug in my code - fixed it. Plots look better now.
#  Learned how to save ggplots, and got the parameters from BMC. Did a ton of
#  tweaking to the plot aesthetics. Hopefully close to done!
#
###############################################################################
## Setup packages
###############################################################################
rm(list=ls()); gc()
startTime <- Sys.time()
set.seed(1029384756)

# Check that required packages are installed
# This does not check that it's the proper version of V2!!!!
if(!all(c('MGDrivE','MGDrivE2','ggplot2') %in% installed.packages()[ ,'Package']) ){
  stop("Packages 'MGDrivE', 'MGDrivE2', and 'ggplot2' are required.")
}


###############################################################################
## Inputs
###############################################################################
####################
## Output Path
####################
baseOut<-'~/Desktop/OUTPUT/BMC'

####################
## Simulation Params
####################
numRep <- 10
tmax <- 3*350
dt <- 1
dt_stoch <- 0.2

relTimes <- seq(from = 100, length.out = 26, by = 7)
relSize <- 100

####################
## Bio and Epi Params
####################
theta <- list(
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
  muF = 0.09,
  muM = 0.09,
  beta = 16,
  nu = 1/(4/24),
  # epidemiological parameters
  NH = 1000,
  X = 0.25,
  f = 1/3,
  Q = 0.9,
  b = 0.05,
  c = 0.15,
  r = 1/200,
  muH = 1/(62*365),
  qEIP = 1/11,
  nEIP = 6
)


###############################################################################
## Inhritance Patterns
###############################################################################
####################
## Cubes
####################
# Simple gene drive
#  Assume equal male/female parameters
#  for simplicity, perfect homing, no female deposition



# cube <- MGDrivE::cubeHomingDrive(cM = 0.5,chM = 1.0,
#                                  cF = 0.5,chF = 1.0,
#                                  dF = 0)
# cube_DIM <- MGDrivE::cubeHomingDrive(cM = 0.5,chM = 1.0,
#                                      cF = 0.5,chF = 1.0,
#                                      dF = 0)


cube <- MGDrivE::cubeMendelian()
cube_DIM <- MGDrivE::cubeMendelian()

####################
## Epi Augmentation
####################
# mosquito to human transmission
m2h <- setNames(object = rep(x = theta$b, times = cube$genotypesN), nm = cube$genotypesID)

# human to mosquito transmission
h2m <- setNames(object = rep(x = theta$c, times = cube$genotypesN), nm = cube$genotypesID)

# add to cubes
cube$b <- m2h
cube$c <- h2m

cube_DIM$b <- m2h
cube_DIM$c <- h2m

####################
## DIM Augmentation
####################
# grab default omega from regular cube
omegaDefault <- cube$omega

# generate new omega
#  pull the original vector
#  for every copy of the gene drive, increase death rate 5x
#   Aa = 5x
#   aa = 10x



# omegaInfect <- cube$omega
# omegaInfect[c("WH","HR","HB")] <- 0.5
# omegaInfect["HH"] <- 0.25


omegaInfect <- cube$omega
omegaInfect["Aa"] <- 0.75
omegaInfect["aa"] <- 0.5





# replace omega from cube with omega for infection
# set in SEI list for infection stages
cube_DIM$omega <- list("S" = omegaDefault,
                       "E" = omegaDefault,
                       "I" = omegaInfect)

####################
## Release of modified mosquitoes
####################
# assume we have a homogenous population in captivity, that we bred for releases
# females
relF <- data.frame("var" = paste0("F_", cube$releaseType, "_", cube$releaseType, "_S"),
                   "time" = relTimes,
                   "value" = relSize,
                   "method" = "add",
                   stringsAsFactors = FALSE)

# males
relM <- data.frame("var" = paste0("M_", cube$releaseType),
                     "time" = relTimes,
                     "value" = relSize,
                     "method" = "add",
                     stringsAsFactors = FALSE)

# combine
releases <- rbind(relF,relM)


###############################################################################
## Build Simulation
###############################################################################
####################
## Set up Petri Net
####################
# Places and transitions
SPN_P <- MGDrivE2::spn_P_epiSIS_node(params = theta, cube = cube)
SPN_T <- MGDrivE2::spn_T_epiSIS_node(spn_P = SPN_P, params = theta, cube = cube)

# Stoichiometry matrix
S <- MGDrivE2::spn_S(spn_P = SPN_P, spn_T = SPN_T)

####################
## Equilibrium
####################
initialCons <- MGDrivE2::equilibrium_SEI_SIS(params = theta, phi = 0.5, log_dd = TRUE,
                                             spn_P = SPN_P, cube = cube)

####################
## Hazards
####################
# using Tau leaping, which is discrete state-space, so exact hazards only
# no DIM
hazards <- MGDrivE2::spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube,
                                 params = initialCons$params, type = "SIS",
                                 log_dd = TRUE, exact = TRUE, verbose = FALSE)

# DIM
hazards_DIM <- MGDrivE2::spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube_DIM,
                                     params = initialCons$params, type = "SIS",
                                     log_dd = TRUE, exact = TRUE, verbose = FALSE)

###############################################################################
## Run Simulation
###############################################################################
# no DIM
PTS_out <- MGDrivE2::sim_trajectory_R(x0 = initialCons$M0, tmax = tmax, dt = dt,
                            dt_stoch = dt_stoch, num_reps = numRep, S = S, hazards = hazards,
                            sampler = "tau", events = releases, verbose = TRUE)

# DIM
PTS_out_DIM <- MGDrivE2::sim_trajectory_R(x0 = initialCons$M0, tmax = tmax, dt = dt,
                                dt_stoch = dt_stoch, num_reps = numRep, S = S, hazards = hazards_DIM,
                                sampler = "tau", events = releases, verbose = TRUE)


###############################################################################
## Analysis and Plotting
###############################################################################
####################
## Analyze
####################
# https://www.statology.org/r-aggregate-multiple-columns/
# https://stackoverflow.com/questions/34523679/aggregate-multiple-columns-at-once
# summarize females/humans by genotype and epi status
fNoDIM <- MGDrivE2::summarize_females_epi(out = PTS_out$state, spn_P = SPN_P)
hNoDIM <- MGDrivE2::summarize_humans_epiSIS(out = PTS_out$state)
hNoDIM$genotype <- "second"
hNoDIM$dim <- "D.I.M. -"

fDIM <- MGDrivE2::summarize_females_epi(out = PTS_out_DIM$state, spn_P = SPN_P)
hDIM <- MGDrivE2::summarize_humans_epiSIS(out = PTS_out_DIM$state)
hDIM$genotype <- "second"
hDIM$dim <- "D.I.M. +"

# aggregate out mosquito genotypes
fNoDIMSEI <- aggregate(x = value ~ time + rep + inf, data = fNoDIM, FUN = sum)
fNoDIMSEI$genotype <- "first"
fNoDIMSEI$dim <- "D.I.M. -"

fDIMSEI <- aggregate(x = value ~ time + rep + inf, data = fDIM, FUN = sum)
fDIMSEI$genotype <- "first"
fDIMSEI$dim <- "D.I.M. +"


####################
## Plot
####################
# data renaming scheme
mNames <- setNames(object = c("Mosquito", "Human"), nm = c("first","second"))

# "Tags" - tag individual plots
labDF <- data.frame(genotype = c("first","first","second","second"),
                   dim = c("D.I.M. -","D.I.M. +","D.I.M. -","D.I.M. +"),
                   label = c("A","B","C","D"))

# y location for the tags
#  Just a note, this is hand-tuned for looks. No programatic way to set it up (that I know of)
yDim <- rep(x = c(max(fNoDIMSEI$value, fDIMSEI$value) * 1.05, max(hNoDIM$value, hDIM$value) * 1.04), each = 2)

# things I use multiple times
axisT <- 18
axisL <- 10
facetL <- 14



"NEED TO ADD TONS OF NOTES TO THIS. I DIDN'T KEEP ANY OF THE SITES I LEARNED FROM."

# plot
curPlot <- ggplot(data = rbind(fNoDIMSEI, fDIMSEI, hNoDIM, hDIM)) +
  geom_line(aes(x = time, y = value, color = inf, group = interaction(rep, inf)), lwd = 0.15) +
  facet_grid(genotype ~ dim, scales = "free_y", labeller = labeller(genotype = mNames)) +
  scale_x_continuous(limits = c(0, tmax), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0,0.14,0.14)) +
  labs(title = waiver(), y = "Count", x = "Time (Days)", color = "Infection \n  Status") +
  geom_text(x = 0.05*tmax, y = yDim, label = labDF$label, data = labDF, fontface = "bold", size = 6) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  theme(axis.title = element_text(face="bold", size = axisT),
        axis.text = element_text(size = axisL),
        strip.text = element_text(face="bold",size = facetL),
        legend.title = element_text(face = "bold", size = facetL),
        legend.text = element_text(size = axisL, hjust = 0),
        legend.key.width = unit(x = 20, units = "pt"),
        plot.margin = unit(c(1,0,-2,1), "pt"),
        legend.margin = margin(t=0, r=0, b=0, l=-10, unit='pt'))

# save output
width <- 170 # full-page width for BMC
aR <- 9/16 # standard wide-screen aspect ratio
ggsave(path = baseOut, filename = "figure_2.tiff",
       plot = curPlot, device = "tiff", compression = "lzw",
       width = width, height = width*aR, units = "mm", dpi = 600)

ggsave(path = baseOut, filename = "figure_2.pdf",
       plot = curPlot, device = "pdf",
       width = width, height = width*aR, units = "mm")

# print time it took
print(paste0('All Sims: ', capture.output(difftime(time1 = Sys.time(), time2 = startTime))))

