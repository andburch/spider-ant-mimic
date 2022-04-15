#devtools::install_github("JimMcL/trajr")
library(trajr)
source("R/functions.R")
source("R/report_results.R")
source("R/tests.R")

# Pre-defined list of categories, in the order they should be reported
allcategories <- c("bee-mimic H. pahangensis", "bee-mimic A. argentifasciata", "bee-mimic P. cruentata",
                   "black stingless bee", "red stingless bee", "dwarf honey bee", "wasp-mimic Pyrophleps", "wasp")

# Any slower than this speed is considered to be hovering
MAX_HOVER_SPEED <- 0.1

# Sizes of plots for manuscript
FIGURE_WIDTH <- 6.267/2*300
FIGURE_HEIGHT <- FIGURE_WIDTH*0.7
FIGURE_RESOLUTION <- 70

# List of info describing calculated parameters
PARAMS_INFO <- list(
  longest_hover_time      = list(in_pca = FALSE, transform = "",     label = "Longest Hover Time (s)", title = "Longest Hover Time"),
  mean_speed              = list(in_pca = TRUE,  transform = "log",  label = "Mean speed (m/s)", title = "Mean speed"),
  min_speed               = list(in_pca = TRUE,  transform = "log",  label = "Minimum speed (m/s)", title = "Minimum speed"),
  max_speed               = list(in_pca = TRUE,  transform = "sqrt", label = "Maximum speed (m/s)", title = "Maximum speed"),
  sd_speed                = list(in_pca = TRUE,  transform = "log",  label = "Standard deviation of speed (m/s)", title = "Standard deviation of speed"),
  sinuosity               = list(in_pca = TRUE,  transform = "log",  label = "Sinuosity", title = "Sinuosity"),
  straightness            = list(in_pca = TRUE,  transform = "",     label = "Straightness", title = "Straightness"),
  first_min_deltaS        = list(in_pca = TRUE,  transform = "",     label = "First min deltaS", title = "First min deltaS"),
  first_min_C             = list(in_pca = TRUE,  transform = "",     label = "First min C", title = "First min C"),
  Emax                    = list(in_pca = TRUE,  transform = "log",  label = "Emax", title = "Emax"),
  directional_change_mean = list(in_pca = TRUE,  transform = "log",  label = "Mean directional change (°/s)", title = "Mean directional change"),
  directional_change_sd   = list(in_pca = TRUE,  transform = "log",  label = "Standard deviation of directional change (°/s)", title = "Standard deviation of directional change")
)

### Functions

longest_hover_time <- function(trj) {
  intervals <- TrajSpeedIntervals(trj, slowerThan = MAX_HOVER_SPEED)
  # Return time of longest interval, or 0 if there are no intervals
  max(c(0, intervals$duration))
}

# Function to calculate the various parameters for a single trajectory
traj_parameters <- function(trj){
  derivs <- TrajDerivatives(trj)
  max_hover <- longest_hover_time(trj)
  mean_speed <- mean(derivs$speed)
  min_speed <- min(derivs$speed)
  max_speed <- max(derivs$speed)
  sd_speed <- sd(derivs$speed)
  sinuosity <- TrajSinuosity(trj)
  straightness <- TrajStraightness(trj)
  resampled <- TrajRediscretize(trj, .001)
  corr <- TrajDirectionAutocorrelations(resampled, round(nrow(resampled) / 4))
  first_min <- TrajDAFindFirstMinimum(corr)
  Emax <- TrajEmax(trj)
  directional_change_mean <- mean(TrajDirectionalChange(trj))
  directional_change_sd <- sd(TrajDirectionalChange(trj))
  # Combine all parameters into a list, and return it
  list(longest_hover_time = max_hover, mean_speed = mean_speed, min_speed = min_speed, max_speed = max_speed, sd_speed = sd_speed, sinuosity = sinuosity, straightness = straightness, first_min_deltaS = first_min[1], first_min_C = first_min[2], Emax = Emax, directional_change_mean = directional_change_mean, directional_change_sd = directional_change_sd)
}

# Transform the various variables so that they are approximately normally distributed.
# The transforms to apply are defined in params_info
normalise_params <- function(params, params_info) {
  norm <- params
  
  for (pn in names(params_info)) {
    pi <- params_info[[pn]]
    if (pi$transform != "")
      # Call the function given its name
      norm[, pn] <- do.call(get(pi$transform), list(norm[, pn]))
  }
  norm
}

#### 
# Read in trajectories from CSV files

# Read in the CSV file which identifies trajectory files together with parameters describing the trajectories
file_list <- read.csv("files_for_analyses.csv", sep = ";", stringsAsFactors = FALSE)
# Omit blank lines
file_list <- na.omit(file_list)
# Read in paths from CSV (actually tab-separated) files, convert to Trajectory objects, scale and smooth them
trjs <- TrajsBuild(file_list$csv.file.name, file_list$fps, file_list$scale, "m", rootDir = "x_y_raw_data", smoothP = 3, smoothN = 101)

# === Calculate trajectory parameters
step_lengths <- TrajsStepLengths(trjs)
mean_step_length <- mean(step_lengths)
params_step1 <- TrajsMergeStats(trjs, traj_parameters)
params <- cbind(file_list, params_step1)

# Build parameters excluding the anomolous wasp mimic which landed part-way through the video
params_excluding_waspmimic <- params[-which(params$csv.file.name == "Pyrophleps_C0073xypts.csv"), ]

# Transform some parameters so they have roughly normal distributions
normalised <- normalise_params(params_excluding_waspmimic, PARAMS_INFO)

# Get a list of species, in the same order as allcategories
all_species <- sapply(allcategories, function(cat) {params$species[which(params$category == cat)[1]]})

# Ensure output directories exist
if (!dir.exists("plots"))
  dir.create("plots")
if (!dir.exists("data_output"))
  dir.create("data_output")


# ==== Output plots and test results

# Plot a selection of trajectories to files
files_to_plot <- c("H.pahangensis_096xypts.csv",
                   "H.pahangensis_C0025xypts.csv",
                   "H.pahangensis_C0041xypts.csv",
                   "H.pahangensis_C0048xypts.csv",
                   "A.argentifasciata_C0139xypts.csv",
                   "Pyrophleps_C0071xypts.csv",
                   "P.cruentata_C0060xypts.csv",
                   "bee_C0016xypts.csv",
                   "red_bee_C0189xypts.csv",
                   "wasp_C0031xypts.csv",
                   "A.andreniformis_C0038xypts.csv")
invisible(sapply(files_to_plot, function(file){
  filename <- file.path("plots", sprintf("%s_trajectory.png", tools::file_path_sans_ext(file)))
  index <- which(params$csv.file.name == file)
  PlotToPng(filename, function() plot(trjs[[index]], main = italicisePattern(params$species[index], ".*"), 
                                      cex.main = 2, cex.axis = 1.3, cex.lab = 2))
}))

# Save lots of results for all trajectories
#report_result(params, "full")

# Save lots of results for all trajectories excluding the anomolous wasp mimic
#report_result(params_excluding_waspmimic, "waspless", dunnTestGroups = params_excluding_waspmimic$species, dunnTestOrder = all_species)

# Report results for normalised variables
report_result(normalised, "normalised")

# Report Kruskall-wallace for (untransformed) parameters not in PCA
for(name in names(PARAMS_INFO)) {
  if (!PARAMS_INFO[[name]]$in_pca)
    report_dunn(params_excluding_waspmimic, name, PARAMS_INFO, "waspless", p.adjustment = "holm")
}

# Produce box plots for selected parameters (un-transformed)
for (name in c("longest_hover_time", "mean_speed", "sinuosity")) {
  PlotToPng(sprintf("plots/bar_plot_%s.png", name), 
            pretty_box_plot(params_excluding_waspmimic[, name], params_excluding_waspmimic$species, all_species, 
                            PARAMS_INFO[[name]]$label, PARAMS_INFO[[name]]$title),
            width = FIGURE_WIDTH, height = FIGURE_HEIGHT, res = FIGURE_RESOLUTION)
}

# Write CSVs with (untransformed) parameter values and summaries
report_params(params_excluding_waspmimic, "waspless", all_species)
