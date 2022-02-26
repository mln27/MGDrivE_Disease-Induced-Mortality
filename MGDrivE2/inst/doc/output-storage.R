## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=7.2, fig.height=4)
set.seed(10)

## -----------------------------------------------------------------------------
# simulation functions
library(MGDrivE2)
# inheritance patterns
library(MGDrivE)
# sparse migration
library(Matrix)

# basic inheritance pattern
cube <- MGDrivE::cubeMendelian()

## -----------------------------------------------------------------------------
# entomological and epidemiological parameters
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
  b = 0.55,
  c = 0.15,
  r = 1/200,
  muH = 1/(62*365),
  qEIP = 1/11,
  nEIP = 2
)

# simulation parameters
tmax <- 50
dt <- 2

## -----------------------------------------------------------------------------
# augment the cube with RM transmission parameters
cube$c <- setNames(object = rep(x = theta$c, times = cube$genotypesN), nm = cube$genotypesID)
cube$b <- c("AA" = theta$b, "Aa" = 0.35, "aa" = 0)

## -----------------------------------------------------------------------------
# nodetypes
node_list <- c("m", "b", "h")
num_nodes <- length(node_list)

# human movement
h_move <- matrix(data = FALSE, nrow = num_nodes, ncol = num_nodes,
                 dimnames = list(node_list, node_list))
h_move[2,3] <- TRUE
h_move[3,2] <- TRUE

# mosquito movement
m_move <- matrix(data = FALSE, nrow = num_nodes, ncol = num_nodes,
                 dimnames = list(node_list, node_list))
m_move[1,2] <- TRUE
m_move[2,1] <- TRUE

# Places and transitions
SPN_P <- spn_P_epiSIS_network(node_list = node_list, params = theta, cube = cube)
SPN_T <- spn_T_epiSIS_network(node_list = node_list, spn_P = SPN_P, params = theta,
                              cube = cube, h_move = h_move, m_move = m_move)

# Stoichiometry matrix
S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)

## -----------------------------------------------------------------------------
# SEI mosquitoes and SIS humans equilibrium
#  outputs required parameters in the named list "params"
#  outputs intial equilibrium for adv users, "init
#  outputs properly filled initial markings, "M0"
initialCons <- equilibrium_SEI_SIS(params = theta, node_list = node_list,
                                   NF = 500, phi = 0.5, NH = 1000, pop_ratio_H = 0.15,
                                   log_dd = TRUE, spn_P = SPN_P, cube = cube)

## -----------------------------------------------------------------------------
# calculate movement rates and movement probabilities
gam <- calc_move_rate(mu = initialCons$params$muF, P = 0.05)

# set mosquito movement rates/probabilities
#  mosquitoes exist in nodes 1 and 2, not 3
mr_mosy <- c(gam, gam, NaN)
mp_mosy <- Matrix::sparseMatrix(i = c(1,2), j = c(2,1), x = 1, dims = dim(m_move))

# set human movement rates/probabilities
#  humans exist in nodes 2 and 3, not 1
mr_human <- c(NaN, 1/7, 1/7)
mp_human <- Matrix::sparseMatrix(i = c(2,3), j = c(3,2), x = 1, dims = dim(h_move))

# put rates and probs into the parameter list
initialCons$params$mosquito_move_rates <- mr_mosy
initialCons$params$mosquito_move_probs <- mp_mosy
initialCons$params$human_move_rates <- mr_human
initialCons$params$human_move_probs <- mp_human

## -----------------------------------------------------------------------------
# exact hazards for integer-valued state space
exact_hazards <- spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube,
                             params = initialCons$params, type = "SIS",
                             log_dd = TRUE, exact = TRUE, tol = NaN,
                             verbose = FALSE)

## -----------------------------------------------------------------------------
# releases
r_times <- seq(from = 10, length.out = 2, by = 7)
r_size <- 50
events <- data.frame("var" = paste0("F_", cube$releaseType, "_", cube$wildType, "_S_1"),
                     "time" = r_times,
                     "value" = r_size,
                     "method" = "add",
                     stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
# main output folder
main_out <- tempdir()

# folders for each stage of analysis
analysis_out <- c("raw", "traces", "analyzed")

# repetitions, 3 of them
rep_out <- formatC(x = 1:3, width=3, format='d', flag='0')

# build analysis folders
analysis_folders <- file.path(main_out, analysis_out)
for(i in analysis_folders){ dir.create(path = i, recursive = TRUE) }

# build repetition folders
rep_folders <- file.path(analysis_folders[1], rep_out)
for(i in rep_folders){ dir.create(i) }

## -----------------------------------------------------------------------------
# delta t
dt_stoch <- 0.1

# run tau-leaping simulation
sim_trajectory_CSV(x0 = initialCons$M0, tmax = tmax, dt = dt,
                   dt_stoch = dt_stoch, folders = rep_folders,
                   stage = c("E", "L", "P", "M", "U", "F", "H"), S = S,
                   hazards = exact_hazards, events = events, verbose = FALSE)

## -----------------------------------------------------------------------------
# list all files from the first repetition
list.files(path = rep_folders[1], full.names = FALSE)

## -----------------------------------------------------------------------------
# read the eggs from repetition 1
egg_stage <- read.csv(file = file.path(rep_folders[1], "E.csv"), header = TRUE)

# read susceptible females from repetition 1
fs_stage <- read.csv(file = file.path(rep_folders[1], "FS.csv"), header = TRUE)

## -----------------------------------------------------------------------------
# eggs
dim(egg_stage)

# susceptible females
dim(fs_stage)

## -----------------------------------------------------------------------------
colnames(egg_stage)

## -----------------------------------------------------------------------------
colnames(fs_stage)

## -----------------------------------------------------------------------------
# split everything by patch, aggregate by genotype
split_aggregate_CSV(read_dir = analysis_folders[1], write_dir = analysis_folders[2],
                    spn_P = SPN_P, tmax = tmax, dt = dt, verbose = FALSE)

## -----------------------------------------------------------------------------
# don't list parent directory
list.dirs(path = analysis_folders[2], recursive = FALSE)

## -----------------------------------------------------------------------------
list.files(path = file.path(analysis_folders[2], rep_out[1]))

## -----------------------------------------------------------------------------
# read the eggs from repetition 1
egg_stage <- read.csv(file = file.path(analysis_folders[2], rep_out[1], "E_0001.csv"),
                      header = TRUE)

# read susceptible females from repetition 1
fs_stage <- read.csv(file = file.path(analysis_folders[2], rep_out[1], "FS_0001.csv"),
                     header = TRUE)

## -----------------------------------------------------------------------------
# eggs
dim(egg_stage)

# susceptible females
dim(fs_stage)

## -----------------------------------------------------------------------------
# eggs
colnames(egg_stage)

# females
colnames(fs_stage)

## -----------------------------------------------------------------------------
# mean and 95% quantiles
summarize_stats_CSV(read_dir = analysis_folders[2], write_dir = analysis_folders[3],
                    spn_P = SPN_P, tmax = tmax, dt = dt, mean = TRUE,
                    quantiles = c(0.025, 0.975), verbose = FALSE)

## -----------------------------------------------------------------------------
list.files(path = analysis_folders[3])

## -----------------------------------------------------------------------------
# read the eggs from repetition 1
egg_stage <- read.csv(file = file.path(analysis_folders[3], "E_Mean_0001.csv"),
                      header = TRUE)

# read susceptible females from repetition 1
fs_stage <- read.csv(file = file.path(analysis_folders[3], "FS_Mean_0001.csv"),
                     header = TRUE)

# eggs
dim(egg_stage)

# susceptible females
dim(fs_stage)

## -----------------------------------------------------------------------------
# eggs
colnames(egg_stage)

# females
colnames(fs_stage)

