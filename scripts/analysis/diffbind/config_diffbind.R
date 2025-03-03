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
  
# 1. DEFINE FILES AND PATHS: ---------------------------------------------------
project <- "atacseq_e16"
analysis <- "diffbind"

input_file <- "data/meta_data/SampleSheet.csv"

output_dir_robj <- file.path("results", "r_objects", analysis)
output_dir_tables <- file.path("results", "tables", analysis)
output_dir_figures <- file.path("results", "figures", analysis)

# Make sure output directory exists, and if it doesn't, create one
if (!dir.exists(output_dir_robj)) {
  dir.create(output_dir_robj, recursive = TRUE)
}

if (!dir.exists(output_dir_tables)) {
  dir.create(output_dir_tables, recursive = TRUE)
}

if (!dir.exists(output_dir_figures)) {
  dir.create(output_dir_figures, recursive = TRUE)
}



# 2. TREATMENT GROUPS: --------------------------------------

# Define control and mutant groups
control_group <- "CTRL_e16"     
mutant_group <- "NICD_e16"  
