# setup_environment.R

# REFERENCE FILE STRUCTURE: ----------------------------------------------------
#project_dir/
#  ├── data/
#  │   ├── meta_data/
#  │   ├── processed_data/
#  ├── logs/
#  ├── notebooks/
#  │   ├── exploratory_analysis/
#  │   ├── final_reports/
#  ├── results/
#  │   ├── figures/
#  │   ├── r_objects/
#  │   ├── reports/
#  │   ├── tables/
#  ├── scripts/
#  │   ├── analysis/
#  │   ├── preprocessing/
#  │   ├── setup/

# ------------------------------------------------------------------------------
# Load renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# Restore the environment from the renv.lock file
if (file.exists("renv.lock")) {
  cat("\nRestoring environment from renv.lock...\n")
  renv::restore()
} else {
  stop("renv.lock file not found. Check that it is in the project directory.")
}

# Set up directory structure
required_directories <- c(
  "data",
  "data/meta_data",
  "data/processed_data",
  "data/motif_files",
  "data/processed_data/peak_files",
  "data/processed_data/bam_files",
  "notebooks",
  "notebooks/exploratory_analysis",
  "notebooks/final_reports",
  "results",
  "results/figures",
  "results/tables",
  "results/r_objects",
  "results/reports",
  "logs",
  "scripts",
  "scripts/analysis",
  "scripts/preprocessing",
  "scripts/setup"
)

for (dir in required_directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("Created directory: ", dir, "\n"))
  }
}

cat("\nEnvironment setup complete.\n")