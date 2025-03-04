# Run after running 'analysis/02_diffbind_e16.R' and 'tables/scripts/tableS3.R'

library(readr)
library(tidyverse)
library(dplyr)

# READ IN DIFFBIND RESULTS: ----------------------------------------------------
diffbind_res_df <- read.csv(
  file.path(output_dir_tables, 
            paste0(project, "_diffbind_results.csv")))

# Create dataframes for CTRL (Fold < 0) keeping only significant peaks 
# (FDR < 0.05)
assign(paste0(control_group, "_enrich"), diffbind_res_df %>% 
         dplyr::filter(FDR < 0.05 & Fold < 0) %>% 
         dplyr::select(seqnames, start, end, peak_id))


# Create dataframes for NICD keeping only significant peaks (FDR < 0.05)
assign(paste0(mutant_group, "_enrich"), diffbind_res_df %>% 
         dplyr::filter(FDR < 0.05 & Fold > 0) %>% 
         dplyr::select(seqnames, start, end, peak_id))


# Write to bed files
filename <- paste0(control_group, "_enrich.bed")

write.table(
  get(paste0(control_group, "_enrich")),   # specify dataframe, not string
  file = file.path(output_dir_tables, filename),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

filename <- paste0(mutant_group, "_enrich.bed")

write.table(
  get(paste0(mutant_group, "_enrich")), # specify dataframe, not string
  file=file.path(output_dir_tables, filename), 
  sep="\t", quote=F, row.names=F, col.names=F)

print(paste0("Bed files saved in ", output_dir_tables))
