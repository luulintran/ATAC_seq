# SET UP

library(tidyverse)
library(dplyr)
library(readr)
library(knitr)

# Set root directory
#knitr::opts_knit$set(root.dir = "/Users/luulitran/Desktop/ATAC_seq/atacseq_nicd_e13_2024_09_18")

# READ IN DATA:------------------------------------------------------------------------------------------
## Read in Diffbind data
CTRL_enriched_list <- read.csv("results/tables/diffbind_tables/CTRL_e13_enriched_list.csv")
NICD_enriched_list <- read.csv("results/tables/diffbind_tables/NICD_e13_enriched_list.csv")

# Read in HOMER findMotifs files for Sox10 (based on de novo motif enrichment) and Sox9 (motif file downloaded from HOMER)
sox10_CTRL <- read.csv("results/tables/homer_tables/SOX10_CTRL_findMotifs_inst.txt", sep ="\t")
sox10_NICD <- read.csv("results/tables/homer_tables/SOX10_NICD_findMotifs_inst.txt", sep ="\t")

sox9_CTRL <- read.csv("results/tables/homer_tables/SOX9_CTRL_findMotifs_inst.txt", sep ="\t")
sox9_NICD <- read.csv("results/tables/homer_tables/SOX9_NICD_findMotifs_inst.txt", sep ="\t")

# CLEAN UP DATAFRAMES: ------------------------------------------------------------------------------------------
# Function to remove unnecessary columns in diffbind dataframes *****************************************************************
process_columns <- function(df) {
  # Remove unnecessary columns
  df <- df[, -c(1, 5:9, 15:19, 21, 22)]
  
  # Change order of columns
  df <- df[, c(7:12, 1:6)]
  
  # Return modified dataframe
  return(df)
}

# Function to clean up homer dataframes **********************************************************************************
process_dataframe <- function(df, tf_name) {
  # Rename the first column to "peak_id"
  colnames(df)[1] <- "peak_id"
  
  # Replace all values in the "Motif.Name" column with the transcription factor name
  df$Motif.Name <- tf_name
  
  # Return cleaned up dataframe
  return(df)
}

# Function to merge dataframes *************************************************************************************************
merge_dataframes <- function(df1, df2) {
  # Merge the data frames based on the "peak_id" column
  merged_df <- merge(df1, df2, by = "peak_id", all = TRUE)
  
  # Remove rows where the "Motif.Name" column is NA or empty
  merged_df <- merged_df[!is.na(merged_df$Motif.Name) & merged_df$Motif.Name != "", ]
  
  # Group by peak_id, keeping duplicates together
  grouped_df <- merged_df %>%
    group_by(peak_id)
  
  # Return the grouped dataframe
  return(grouped_df)
}
# ********************************************************************************************************************************
## Clean up dataframes
## Apply the function to remove unnecessary columns and reorder the remaining columns
CTRL_enriched_list <- process_columns(CTRL_enriched_list)
NICD_enriched_list <- process_columns(NICD_enriched_list)

## Apply function to replace column name of PositionID and values in Motif.Name
sox10_CTRL <- process_dataframe(sox10_CTRL, "Sox10")
sox9_CTRL <- process_dataframe(sox9_CTRL, "Sox9")
sox10_NICD <- process_dataframe(sox10_NICD, "Sox10")
sox9_NICD <- process_dataframe(sox9_NICD, "Sox9")

## Merge dataframes into one, remove NA's (peaks that didn't have the motif), and group duplicate peaks (peaks that contain the motif more than once, depending on position)
# Apply merge function. Now you should have a dataframe containing the peaks that contain the motif
sox10_CTRL_motif_inst <- merge_dataframes(CTRL_enriched_list, sox10_CTRL)
sox9_CTRL_motif_inst <- merge_dataframes(CTRL_enriched_list, sox9_CTRL)
sox10_NICD_motif_inst <- merge_dataframes(NICD_enriched_list, sox10_NICD)
sox9_NICD_motif_inst <- merge_dataframes(NICD_enriched_list, sox9_NICD)

# Reorder columns
sox10_CTRL_motif_inst <- sox10_CTRL_motif_inst[, c(1, 7:9, 2:6, 10:ncol(sox10_CTRL_motif_inst))]
sox9_CTRL_motif_inst <- sox9_CTRL_motif_inst[, c(1, 7:9, 2:6, 10:ncol(sox9_CTRL_motif_inst))]
sox10_NICD_motif_inst <- sox10_NICD_motif_inst[, c(1, 7:9, 2:6, 10:ncol(sox10_NICD_motif_inst))]
sox9_NICD_motif_inst <- sox9_NICD_motif_inst[, c(1, 7:9, 2:6, 10:ncol(sox9_NICD_motif_inst))]

# Reorder by FDR
sox10_CTRL_motif_inst <- sox10_CTRL_motif_inst[order(sox10_CTRL_motif_inst$FDR), ]
sox9_CTRL_motif_inst <- sox9_CTRL_motif_inst[order(sox9_CTRL_motif_inst$FDR), ]

sox10_NICD_motif_inst <- sox10_NICD_motif_inst[order(sox10_NICD_motif_inst$FDR), ]
sox9_NICD_motif_inst <- sox9_NICD_motif_inst[order(sox9_NICD_motif_inst$FDR), ]

# WRITE TO FILES: ----------------------------------------------------------------------------------------------------
write.csv(sox10_CTRL_motif_inst, "results/tables/homer_tables/sox10_CTRL_motif_inst_annotated.csv")
write.csv(sox9_CTRL_motif_inst, "results/tables/homer_tables/sox9_CTRL_motif_inst_annotated.csv")

write.csv(sox10_NICD_motif_inst, "results/tables/homer_tables/sox10_NICD_motif_inst_annotated.csv")
write.csv(sox9_NICD_motif_inst, "results/tables/homer_tables/sox9_NICD_motif_inst_annotated.csv")

print("HOMER motifs instances merged with Diffbind results. Results are saved in the 'results/tables/' directory.")

