# Collects all trajectory CSV files into a directory
OUTDIR <- "trajectories"
INDIR <- "x_y_raw_data" 
file_list <- read.csv("files_for_analyses.csv", sep = ";", stringsAsFactors = FALSE)
file_list <- na.omit(file_list)

if (!dir.exists(OUTDIR))
  dir.create(OUTDIR)
for (src in file_list$csv.file.name) {
  file <- list.files(INDIR, pattern = src, full.names=TRUE, recursive=TRUE)
  #cat(sprintf("%s => '%s'\n", src, file))
  file.copy(file, OUTDIR)
}