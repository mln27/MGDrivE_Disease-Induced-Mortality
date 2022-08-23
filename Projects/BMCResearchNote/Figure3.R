###############################################################################
#      ____  __  _________            _____
#     / __ )/  |/  / ____/           |__  /
#    / __  / /|_/ / /      ______     /_ <
#   / /_/ / /  / / /___   /_____/   ___/ /
#  /_____/_/  /_/\____/            /____/
# .............................................................................
# BMC Research Note - Figure 3
# Ndeffo Mbah Lab
# Jared Bennett, jared_bennett@berkeley.edu
#
# 20220527
#  This script is designed to generate figure 3 of the BMC research not.
#  The works stems significantly from the life-cycle node vignette. This is
#  augmented with an SEM-cis cube to demonstrate the new cube, the small-molecule
#  style releases, and the new analysis function.
#  Again, the point is simply to show that the idea and implementation work,
#  NOT to demonstrate any kind of scientific point.
#
#
###############################################################################
## Setup packages
###############################################################################
rm(list=ls()); gc()
startTime <- Sys.time()
set.seed(1029384756)

# Check that required packages are installed
# This does not check that it's the proper version of V2!!!!
if(!all(c('MGDrivE2','ggplot2','parallel','gridExtra') %in% installed.packages()[ ,'Package']) ){
  stop("Packages 'MGDrivE2', 'ggplot2', 'parallel', and 'gridExtra' are required.")
}

library(ggplot2)


###############################################################################
## Inputs
###############################################################################
####################
## Output Path
####################
baseOut <-'~/Desktop/OUTPUT/BMC'
scratch <- file.path(baseOut, 'scratchF3')

####################
## Simulation Params
####################
numCores <- 2
numRep <- 10
tmax <- 2*365
dt <- 1
dt_stoch <- 0.1

repNames <- formatC(x = 1:numRep, width = 3, format = 'd', flag = '0')
runFolders <- file.path(scratch, repNames)

if(!dir.exists(paths = scratch)) dir.create(path = scratch, recursive = FALSE)
for(i in runFolders) dir.create(path = i)


####################
## Releases/Spraying
####################
# gene drive releases
#  start at ~3 months in, release weekly for 1 month
gdRelTimes <- seq(from = 90, length.out = 4, by = 7)
gdRelSize <- 100

# small-molecule spraying to turn on SEM construct
#  give GD ~9mo to spread, then turn on SEM and spray for ~6mo
SEMRelTimes <- seq(from = 365, length.out = 182, by = 1)
numSEMRT <- length(SEMRelTimes)

# Small-molecule sprays are setup to switch a genotype from "off" to "on"
# This requires the user to specify "off" and "on" genotypes for the algorithm.
# These are specific to the inheritance pattern.
#  the first row is the "off" genotypes - just GD, no SEM
#  the second row is what those genotypes become when turned "on" - GE and SEM
#  these have to be in the genotypes of the cube being simulated
swapMat <- matrix(data = c("GW","GG","GU","GR","GV","GH","GS",
                           "HW","HH","HU","HR","HV","HH","HS"),
                  nrow = 7, ncol = 2, dimnames = list(NULL, c("off","on")))

####################
## Bio and Epi Params
####################
# entomological parameters
#  from lifecycle node vignette
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
  nu = 1/(4/24),
  phi = 0.5
)

NF <- 1000


###############################################################################
## Inheritance Patterns
###############################################################################
####################
## Cube
####################
# original, cis-acting SEM construct
#  equal parameters between males and females
#  all deposition is turned off
#  no genotype-specific fitness costs
#  no resistance generation - this is not realistic, but simpler to interpret here
cube <- MGDrivE2::cubeSEMcis(pF = 0.5, qF = 1.0, rF = 0.0,
                             aF = 0.5, bF = 1.0, cF = 1.0)


###############################################################################
## Build Simulation
###############################################################################
####################
## Set up Petri Net
####################
# Places and transitions
SPN_P <- MGDrivE2::spn_P_lifecycle_node(params = theta, cube = cube)
SPN_T <- MGDrivE2::spn_T_lifecycle_node(spn_P = SPN_P, params = theta,
                                        cube = cube, feqTol = 1e-10)

# Stoichiometry matrix
S <- MGDrivE2::spn_S(spn_P = SPN_P, spn_T = SPN_T)

####################
## Equilibrium
####################
initialCons <- MGDrivE2::equilibrium_lifeycle(params = theta, NF = NF,
                                              phi = theta$phi, log_dd = TRUE,
                                              spn_P = SPN_P, cube = cube)

####################
## Hazards
####################
# using Tau leaping, which is discrete state-space, so exact hazards only
hazards <- MGDrivE2::spn_hazards(spn_P = SPN_P, spn_T = SPN_T, cube = cube,
                                 params = initialCons$params, type = 'life',
                                 log = TRUE, exact = TRUE, verbose = FALSE)

####################
## Releases and Spraying
####################
# Male releases only - this is standard in mosquito field because females bite
#  but males do not

# gene drive releases
gdEvents <- data.frame('var' = paste0('M_', cube$releaseType),
                       'time' = gdRelTimes,
                       'value' = gdRelSize,
                       'method' = 'add',
                       stringsAsFactors = FALSE)

# small-molecule spraying
# this is how the spn builds the statespace
#  we have to recreate this perfectly, or the swapping won't work
# females_unmated <- file.path("U",g, fsep = "_")
# females_unmated <- file.path("U_", g, "_", node_id, fsep = "")
#
# males <- file.path("M",g,fsep = "_")
# males <- file.path("M",g,node_id,fsep = "_")
#
# females <- file.path("F_",rep(g, each = nG),"_",g, fsep = "")
# females <- file.path("F_",rep(g, each = nG),"_",g,"_",node_id, fsep = "")
uEvents <- data.frame('var' = rep(paste0("U_", swapMat[ ,"on"]), each = numSEMRT),
                      'value' = rep(paste0("U_", swapMat[ ,"off"]),  each = numSEMRT),
                      'time' = SEMRelTimes,
                      'method' = "swap")

mEvents <- data.frame('var' = rep(paste0("M_", swapMat[ ,"on"]), each = numSEMRT),
                      'value' = rep(paste0("M_", swapMat[ ,"off"]), each = numSEMRT),
                      'time' = SEMRelTimes,
                      'method' = "swap")

# grab all female states from petri net
femHold <- grep(pattern = "^F_", x = SPN_P$u, value = TRUE)
# grab only the states with impacted genotypes
femVal <- grep(pattern = paste0(swapMat[ ,"off"], collapse = "|"), x = femHold, value = TRUE)
# swap appropriate genotypes
femVar <- femVal
for(i in 1:NROW(swapMat)){
  femVar <- gsub(pattern = swapMat[i,"off"], replacement = swapMat[i,"on"], x = femVar, fixed = TRUE)
}
# female events
fEvents <- data.frame('var' = rep(femVar, each = numSEMRT),
                      'value' = rep(femVal, each = numSEMRT),
                      'time' = SEMRelTimes,
                      'method' = "swap")

# combine all events for simulation
releases <- rbind(gdEvents, uEvents, mEvents, fEvents)


###############################################################################
## Run Simulation
###############################################################################
####################
# Sims
####################
# Setup the cluster
#  https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
if(Sys.info()[['sysname']] == "Windows"){
  # windows can't run fork clusters
  #  it has to use sockets
  cl <- parallel::makePSOCKcluster(names = numCores)

  # load libraries on each socket
  #  if we call using package::function(), this is unnecessary
  # parallel::clusterEvalQ(cl = cl, expr = {library(MGDrivE2)})

  # export required objects to each socket
  parallel::clusterExport(cl = cl,
                          varlist = c('initialCons', 'tmax', 'dt',
                                      'dt_stoch', 'runFolders', 'S',
                                      'hazards', 'releases') )

} else {
  # *nix systems
  # Sys.info() does distinguish Linux vs Mac, but we don't care for clusters
  # no copying, no memory issues!
  cl <- parallel::makeForkCluster(nnodes = numCores)

}

# Set parallel seed
parallel::clusterSetRNGStream(cl = cl, iseed = sample(x = (-.Machine$integer.max):.Machine$integer.max,
                                                      size = 1, replace = FALSE))

# run reps
#  "retList" just hides the return value, since the function doesn't return data
retList <- parallel::clusterApplyLB(cl = cl, x = 1:numRep, fun = function(x){
  MGDrivE2::sim_trajectory_CSV(
    x0 = initialCons$M0, tmax = tmax, dt = dt, dt_stoch = dt_stoch,
    folders = runFolders[x], sampler = "tau", stage = c("M", "F"), S = S,
    hazards = hazards, events = releases, verbose = FALSE, maxhaz = 1e12
  )

  gc()
})

# stop cluster
parallel::stopCluster(cl)

####################
# Aggregate
####################
MGDrivE2::split_aggregate_CSV(read_dir = scratch, spn_P = SPN_P,
                              tmax = tmax, dt = dt, stage = c("M","FS"),
                              erlang = FALSE, sum_fem = FALSE,
                              verbose = FALSE, rem_file = TRUE)


###############################################################################
## Analysis and Plotting
###############################################################################
####################
## Analyze
####################
# genotype list
goiList <- c(setNames(object = as.list(cube$genotypesID), nm = cube$genotypesID),
             list("Total"=cube$genotypesID))

# allele list
#  these are all the alleles in this cube
aVec <- c("W","G","H","U","R","V","S")

holdList <- lapply(X = aVec, FUN = function(x){
  c(grep(pattern = paste0(x,'.|.',x), x = cube$genotypesID, value = TRUE),
    grep(pattern = paste0(x,x), x = cube$genotypesID, value = TRUE))
})

# the "Total" is kinda a cheat - I know it's diploid, therefore the total
#  allele count is just 2 * number of genotypes
alleleList <- c(setNames(object = holdList, nm = aVec),
                list("Total" = rep.int(x = cube$genotypesID, times = 2)))

# do analysis
MGDrivE2::analyze_ggplot_CSV(read_dir = scratch, name = "genotypes", sex = "agg",
                             patch_agg = FALSE, goi = goiList, drop_zero_goi = TRUE)

MGDrivE2::analyze_ggplot_CSV(read_dir = scratch, name = "alleles", sex = "agg",
                             patch_agg = FALSE, goi = alleleList, drop_zero_goi = TRUE)


####################
## Plot
####################
# references
# https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/
# https://stackoverflow.com/questions/44170871/how-does-ggplot-scale-continuous-expand-argument-work
# https://ggplot2.tidyverse.org/reference/theme.html
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# read in data
gDat <- read.csv(file = file.path(scratch, "analysis_genotypes/metrics.csv"),
                 header = TRUE, sep = ",", check.names = FALSE)
aDat <- read.csv(file = file.path(scratch, "analysis_alleles/metrics.csv"),
                 header = TRUE, sep = ",", check.names = FALSE)

# things I use multiple times
axisT <- 14
axisL <- 8
facetL <- 10


# data + color scheme and organization
gPlot <- ggplot(data = gDat, aes(x = Time, y = Mean, group = interaction(GOI), color = GOI)) +
  # facet, for label only
  facet_grid(Patch ~ Sex, scales = "free_y", labeller = labeller(Sex = c("agg" = "Genotypes"))) +
  # ribbon and fill
  geom_ribbon(data = gDat, aes(ymin = get("2.5%"), ymax = get("97.5%"), fill = GOI),
              alpha = 0.25, size = 0.0, show.legend = FALSE) +
  # mean
  geom_line(data = gDat, mapping = aes(x = Time, y = Mean,
                                       group = interaction(GOI),
                                       color = GOI), size = 0.3) +
  ## aesthetics
  # color
  scale_color_hue(h = c(-80,80)) +
  scale_fill_hue(h = c(-80,80)) +
  # remove extra space at edges of graph
  scale_x_continuous(limits = c(0, tmax), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0,0.08,0)) +
  # fix labels
  labs(title = waiver(), y = "Genotype Count", x = "Time (Days)") +
  # label plot
  geom_text(x = 0.05*tmax, y = max(gDat$Mean)*1.05, label = "A",
            fontface = "bold", size = 7, inherit.aes = FALSE) +
  # change theme
  theme_bw() +
  theme(
    # remove Y facet label
    strip.text.y = element_blank(),
    # x/y axis
    axis.title = element_text(face="bold", size = axisT),
    axis.text = element_text(size = axisL),
    # facet label
    strip.text = element_text(face="bold",size = facetL),
    # legend
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = facetL),
    legend.text = element_text(size = axisL),
    legend.spacing.x = unit(x = 1, units = "pt"),
    legend.spacing.y = unit(x = 1, units = "pt"),
    legend.margin = margin(t=-10, r=0, b=0, l=0, unit='pt'),
    # plot margins
    plot.margin = unit(c(1,4,0,1), "pt")
  ) +
  # more legend work
  guides(colour = guide_legend(title = "Genotypes", title.position = "top",
                               title.hjust = 0.5, override.aes = list(size=3),
                               keywidth = unit(x = 20, units = "pt"),
                               nrow = 1, label.position = "top",
                               keyheight = unit(x = 1, units = "pt")))


# data + color scheme and organization
aPlot <- ggplot(data = aDat, aes(x = Time, y = Mean, group = interaction(GOI), color = GOI)) +
  # facet, for label only
  facet_grid(Patch ~ Sex, scales = "free_y", labeller = labeller(Sex = c("agg" = "Alleles"))) +
  # ribbon and fill
  geom_ribbon(data = aDat, aes(ymin = get("2.5%"), ymax = get("97.5%"), fill = GOI),
              alpha = 0.25, size = 0.0, show.legend = FALSE) +
  # mean
  geom_line(data = aDat, mapping = aes(x = Time, y = Mean,
                                       group = interaction(GOI),
                                       color = GOI), size = 0.3) +
  ## aesthetics
  # color
  scale_color_hue(h = c(100,260)) +
  scale_fill_hue(h = c(100,260)) +
  # remove extra space at edges of graph
  scale_x_continuous(limits = c(0, tmax), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0,0.08,0)) +
  # fix labels
  labs(title = waiver(), y = "Allele Count", x = "Time (Days)") +
  # label plot
  geom_text(x = 0.05*tmax, y = max(aDat$Mean)*1.05, label = "B",
            fontface = "bold", size = 7, inherit.aes = FALSE) +
  # change theme
  theme_bw() +
  theme(
    # remove Y facet label
    strip.text.y = element_blank(),
    # x/y axis
    axis.title = element_text(face="bold", size = axisT),
    axis.text = element_text(size = axisL),
    # facet label
    strip.text = element_text(face="bold",size = facetL),
    # legend
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = facetL),
    legend.text = element_text(size = axisL),
    legend.spacing.x = unit(x = 5, units = "pt"),
    legend.spacing.y = unit(x = 1, units = "pt"),
    legend.margin = margin(t=-10, r=0, b=0, l=0, unit='pt'),
    # plot margins
    plot.margin = unit(c(1,1,0,4), "pt")
  ) +
  # more legend work
  guides(colour = guide_legend(title = "Alleles", title.position = "top",
                               title.hjust = 0.5, override.aes = list(size=3),
                               keywidth = unit(x = 20, units = "pt"),
                               nrow = 1, label.position = "top",
                               keyheight = unit(x = 1, units = "pt")))


# save output
width <- 170 # full-page width for BMC
aR <- 9/16 # standard wide-screen aspect ratio
ggsave(path = baseOut, filename = "figure_3.tiff",
       plot = gridExtra::grid.arrange(gPlot, aPlot, ncol = 2), device = "tiff", compression = "lzw",
       width = width, height = width*aR, units = "mm", dpi = 600)

ggsave(path = baseOut, filename = "figure_3.pdf",
       plot = gridExtra::grid.arrange(gPlot, aPlot, ncol = 2), device = "pdf",
       width = width, height = width*aR, units = "mm")



####################
## Cleanup
####################
# delete data
unlink(x = scratch, recursive = TRUE)

# print time it took
print(paste0('All Sims: ', capture.output(difftime(time1 = Sys.time(), time2 = startTime))))


