# SET UP
library(DiffBind)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(dplyr)
library(readr)

# LOAD RDS FILE OF DIFFBIND DBA OBJECT: ----------------------------------------
dbObj <- readRDS(
  file = file.path(output_dir_robj, "diffbind_dbObj.rds"))

# Create dba report with results including stats
dbObj.DB <- dba.report(dbObj)

# Make it into a dataframe and save for later
out <- as.data.frame(dbObj.DB)

# Write to file
filename <- "diffbind_deseq2_results.txt"
write.table(out, 
            file=file.path(output_dir_tables, filename), 
            sep="\t", 
            quote=F, 
            row.names=F)

print(paste0(filename, " saved in ", output_dir_tables))