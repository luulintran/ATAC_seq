# SET UP
library(DiffBind)
library(tidyverse)
library(dplyr)
library(readr)

# Set working directory to project
# setwd("..")

# READ IN DATA:-----------------------------------------------------------------
# Sample sheet containing paths to filtered bam files and narrowpeak files

samples <- read.csv(input_file)


# CREATE DBA OBJECT: -----------------------------------------------------------
dbObj <- dba(sampleSheet=samples)

# COUNT READS: -----------------------------------------------------------------
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, summits = 75)


# NORMALIZE DATA: --------------------------------------------------------------
# Normalization is based on sequencing depth by default
dbObj <- dba.normalize(dbObj)

# Save the DBA object after normalizing read counts for later use (making normalized counts plots)
filename <- "norm_read_counts_dbObj.rds"
saveRDS(dbObj, file = file.path(output_dir_robj,filename))

# ESTABLISH MODEL DESIGN AND CONTRAST: -----------------------------------------
# Make the "Control" condition the baseline or denominator of the default contrast.
dbObj <- dba.contrast(dbObj, 
                      reorderMeta = list(Condition = control_group))

# PERFORM DIFFERENTIAL ANALYSIS: -----------------------------------------------
# By default, using DESeq2 but can change to edgeR
dbObj <- dba.analyze(dbObj)

# SAVE DBA OBJECT FOR LATER USE: -----------------------------------------------
# Save the DBA object after differential analysis
filename <- "diffbind_dbObj.rds"
saveRDS(dbObj, file = file.path(output_dir_robj, filename))

print(paste0("Diffbind analysis done. RDS files saved in ", output_dir_robj))


