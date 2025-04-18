
# Run after running 'analysis/02_diffbind_e16.R'

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

# CONVERT SEQNAMES TO UCSC STYLE: ----------------------------------------------
# We want UCSC style and keep mm10 genome build consistent
# Convert to UCSC style from Ensembl

# Takes the seqnames of GRanges object and makes a vector of new seqnames that
# matches the UCSC style.
UCSC_newstyle <- mapSeqlevels(seqlevels(dbObj.DB), "UCSC") 

# Changes the seqnames in dbObj to UCSC style, so now seqnames column has "chr"
dbObj.DB_UCSC <- renameSeqlevels(dbObj.DB, UCSC_newstyle) 

# Store as a dataframe
out <- as.data.frame(dbObj.DB_UCSC)

# Add peak id as "peak_(rownumber)" for later use, 
# like for HOMER motif enrichment analysis
out <- out %>%
  mutate(peak_id = paste("peak_", seqnames, "_",row_number(), sep = ""))

# ASSIGN PEAKS TO NEAREST TSS: -------------------------------------------------
# Using ChIPseeker

# Get annotation data
# This is the UCSC annotation data for mm10
# The ATAC-seq data was mapped to mm10 so we want to keep it consistent:
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate differential peaks to nearest TSS ***********************************
peak_gr <- GRanges(seqnames = out$seqnames,
                   ranges = IRanges(start = out$start, end = out$end))

# Annotate +/- 3 kb around the TSS
peak_anno <- annotatePeak(peak_gr, 
                          tssRegion = c(-3000,3000), 
                          TxDb = txdb, 
                          annoDb = "org.Mm.eg.db")
# Store annotated peaks as a dataframe
peak_anno_df <- as.data.frame(peak_anno)

# merge the out df (containing stats) and the 
# peak_anno_df (containing annotations) and match by seqnames, start, and end
merged_df <- merge(out, peak_anno_df, by = c("seqnames", "start", "end"))
merged_df <- subset(merged_df, select = -c(width.y, strand.y))

# Order df by FDR
merged_df <- merged_df[order(merged_df$FDR), ]

# Write annotated diffbind results to file which can be opened in excel
filename <- paste0(project,"_diffbind_results.csv")
write.csv(
  merged_df, 
  file= file.path(output_dir_tables, filename), 
  row.names=F)

print(paste0(filename, " saved in ", output_dir_tables))

