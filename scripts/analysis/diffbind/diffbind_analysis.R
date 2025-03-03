# SET UP
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = FALSE)
BiocManager::install(c("DiffBind", 
                       "profileplyr", 
                       "ChIPseeker", 
                       "org.Mm.eg.db", 
                       "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                       "pheatmap"), ask = FALSE)

library(DiffBind)
library(profileplyr)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(dplyr)
library(readr)

# Set working directory to project
# setwd("..")

# READ IN DATA:------------------------------------------------------------------------------------------
# Sample sheet containing paths to filtered bam files and narrowpeak files
samples <- read.csv("data/meta/SampleSheet.csv")

# CREATE DBA OBJECT: -------------------------------------------------------------------------------------
dbObj <- dba(sampleSheet=samples)

# COUNT READS: -------------------------------------------------------------------------------------------
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, summits = 75)

# FRiP scores
info <- dba.show(dbObj)
libsizes <- cbind(LibReads=info$Reads, FRiP =info$FRiP,
                  PeakReads = round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes

# Correlation sample heatmap
filename <- "results/figures/diffbind_figures/corr_heatmap.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plot(dbObj)

dev.off()

# NORMALIZE DATA: -------------------------------------------------------------------------------------------
# Normalization is based on sequencing depth by default
dbObj <- dba.normalize(dbObj)

# ESTABLISH MODEL DESIGN AND CONTRAST: ---------------------------------------------------------------------
# Make the "Control" condition the baseline or denominator of the default contrast.
dbObj <- dba.contrast(dbObj, 
                      reorderMeta = list(Condition = "Control"))

# PERFORM DIFFERENTIAL ANALYSIS: ---------------------------------------------------------------------
# By default, using DESeq2 but can change to edgeR
dbObj <- dba.analyze(dbObj)
dba.show(dbObj, bContrasts = TRUE)

# Correlation heatmap of differentially accessible sites
filename <- "results/figures/diffbind_figures/corr_heatmap_diffsites.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plot(dbObj, contrast = 1)

dev.off()

# SAVE DBA OBJECT FOR LATER USE: -------------------------------------------------------------------------
# Save the DBA object after differential analysis
saveRDS(dbObj, file = "results/r_objects/diffbind_dbObj.rds")

# RETRIEVE DIFFBIND RESULTS: ------------------------------------------------------------------------------
dbObj.DB <- dba.report(dbObj)

# EXPLORATORY ANALYSIS: -----------------------------------------------------------------------------------

# Venn Diagram
filename <- "results/figures/diffbind_figures/venn_diagram.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotVenn(dbObj, contrast = 1, bDB = TRUE,
             bGain = TRUE, bLoss = TRUE, bAll = FALSE)

dev.off()

# PCA Plot ******************************************************

# PCA on all sites
filename <- "results/figures/diffbind_figures/pca_plot_allsites.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotPCA(dbObj, DBA_CONDITION, label = DBA_ID)

dev.off()

# PCA on differential sites
filename <- "results/figures/diffbind_figures/pca_plot_diffsites.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotPCA(dbObj, contrast = 1, label = DBA_CONDITION)

dev.off()

# MA plot ******************************************************
filename <- "results/figures/diffbind_figures/MA_plot.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotMA(dbObj)

dev.off()

# Volcano plot ******************************************************
filename <- "results/figures/diffbind_figures/volcano_plot.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotVolcano(dbObj)

dev.off()

# Box plot ******************************************************
filename <- "results/figures/diffbind_figures/box_plot.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotBox(dbObj)

dev.off()

# Heatmap  ******************************************************

# Heatmap colors
filename <- "results/figures/diffbind_figures/heatmap_diffsites.png"
png(filename, width = 5, height = 5, units = "in", res = 300)

hmap <- colorRampPalette(c("white", "#938fc2", "#40027e")) (n =100)
readscores <- dba.plotHeatmap(dbObj, contrast = 1, correlations = FALSE,
                              scale = "row", colScheme = hmap)
print(readscores)

dev.off()

# Profile plot ******************************************************
filename <- "results/figures/diffbind_figures/profile_plot.png"
png(filename, width = 5, height = 5, units = "in", res = 300)

profiles <- dba.plotProfile(dbObj, merge = DBA_REPLICATE)
dba.plotProfile(profiles)

dev.off()

# SAVE DIFFBIND RESULTS AS TXT FILE: -----------------------------------------------------------------------
# Without gene annotation
# This is the GRanges object containing info on the differential sites including stats
dbObj.DB <- dba.report(dbObj)

# Turn it into a dataframe
out <- as.data.frame(dbObj.DB)

# Write to file
write.table(out, file="results/tables/diffbind_tables/diffbind_deseq2_results.txt", sep="\t", quote=F, row.names=F)

# CONVERT SEQNAMES TO UCSC STYLE: ------------------------------------------------------------------------------
#We want UCSC style and keep mm10 genome build consistent
#Convert to UCSC style from Ensembl

#Takes the seqnames of GRanges object and makes a vector of new seqnames that match the UCSC style
UCSC_newstyle <- mapSeqlevels(seqlevels(dbObj.DB), "UCSC") 

#Changes the seqnames in dbObj to UCSC style, so now seqnames column has "chr"
dbObj.DB_UCSC <- renameSeqlevels(dbObj.DB, UCSC_newstyle) 

# MAKE BED FILES OF THE ENRICHED PEAKS: ---------------------------------------------------------------------
out <- as.data.frame(dbObj.DB_UCSC)

# Create dataframes for CTRL keeping only significant peaks (FDR < 0.05)
CTRL_enrich <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end)

#add peak id as "peak_(row number)"
CTRL_enrich <- CTRL_enrich %>%
  mutate(peak_id = paste("peak_", seqnames, "_",row_number(), sep = ""))

# This dataframe has everything including the stats info
CTRL_enrich_all <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0)

CTRL_enrich_all <- CTRL_enrich_all %>%
  mutate(peak_id = paste("peak_", seqnames, "_",row_number(), sep = ""))

# Write to file
write.table(CTRL_enrich, file="results/tables/diffbind_tables/CTRL_e13_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(CTRL_enrich_all, file="results/tables/diffbind_tables/CTRL_e13_enriched_allstats.bed", sep="\t", quote=F, row.names=F, col.names=F)



NICD_enrich <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end)

#add peak id as "peak_(row number)"
NICD_enrich <- NICD_enrich %>%
  mutate(peak_id = paste("peak_", seqnames, "_",row_number(), sep = ""))


NICD_enrich_all <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0)

NICD_enrich_all <- NICD_enrich_all %>%
  mutate(peak_id = paste("peak_", seqnames, "_",row_number(), sep = ""))

# Write to file
write.table(NICD_enrich, file="results/tables/diffbind_tables/NICD_e13_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(NICD_enrich_all, file="results/tables/diffbind_tables/NICD_e13_enriched_allstats.bed", sep="\t", quote=F, row.names=F, col.names=F)

# ASSIGN PEAKS TO NEAREST TSS: ------------------------------------------------------------------------------
# Using ChIPseeker

# Get annotation data
#This is the UCSC annotation data for mm10
# The ATAC-seq data was mapped to mm10 so want to keep it consistent:
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate CTRL peaks to nearest TSS ***********************************************************
CTRL_peak_gr <- GRanges(seqnames = CTRL_enrich_all$seqnames,
                        ranges = IRanges(start = CTRL_enrich_all$start, end = CTRL_enrich_all$end))

# Annotate +/- 3 kb around the TSS
CTRL_peak_anno <- annotatePeak(CTRL_peak_gr, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
CTRL_peak_anno_df <- as.data.frame(CTRL_peak_anno)

# merge dfs and match seqnames, start, and end
CTRL_merged_df <- merge(CTRL_enrich_all, CTRL_peak_anno_df, by = c("seqnames", "start", "end"))
CTRL_merged_df <- subset(CTRL_merged_df, select = -c(width.y, strand.y))

#Order df by FDR
CTRL_merged_df <- CTRL_merged_df[order(CTRL_merged_df$FDR), ]

# Annotate NICD peaks to nearest TSS *************************************************************
# NICD peaks GRanges object
NICD_peak_gr <- GRanges(seqnames = NICD_enrich_all$seqnames,
                        ranges = IRanges(start = NICD_enrich_all$start, end = NICD_enrich_all$end))

# Annotate +/- 3 kb around the TSS
NICD_peak_anno <- annotatePeak(NICD_peak_gr, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
NICD_peak_anno_df <- as.data.frame(NICD_peak_anno)

# merge dfs and match seqnames, start, and end
NICD_merged_df <- merge(NICD_enrich_all, NICD_peak_anno_df, by = c("seqnames", "start", "end"))
NICD_merged_df <- subset(NICD_merged_df, select = -c(width.y, strand.y))

#Order df by FDR
NICD_merged_df <- NICD_merged_df[order(NICD_merged_df$FDR), ]


#rename dataframes
NICD_enriched_list <- NICD_merged_df
CTRL_enriched_list <- CTRL_merged_df

cat("Number of gained peaks:", nrow(NICD_enriched_list), "\n")
cat("Number of lost peaks:", nrow(CTRL_enriched_list), "\n")

# Save enriched peaks as CSV files
write.csv(NICD_enriched_list, "results/tables/diffbind_tables/NICD_e13_enriched_list.csv")
write.csv(CTRL_enriched_list, "results/tables/diffbind_tables/CTRL_e13_enriched_list.csv")


# VISUALIZE GENOMIC FEATURES: ------------------------------------------------------------------------------

# Annotation pie chart *********************************************************
filename <- "results/figures/chipseeker_figures/CTRL_anno_pie.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plotAnnoPie(CTRL_peak_anno)

dev.off()

filename <- "results/figures/chipseeker_figures/NICD_anno_pie.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plotAnnoPie(NICD_peak_anno)

dev.off()

# Annotation bar chart ********************************************************
filename <- "results/figures/chipseeker_figures/CTRL_anno_bar.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plotAnnoBar(CTRL_peak_anno)

dev.off()

filename <- "results/figures/chipseeker_figures/NICD_anno_bar.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plotAnnoBar(NICD_peak_anno)

dev.off()

# Distance to TSS *************************************************************
filename <- "results/figures/chipseeker_figures/CTRL_tss_dist.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plotDistToTSS(CTRL_peak_anno,
              title="Distribution of accessibility\nrelative to TSS")

dev.off()

filename <- "results/figures/chipseeker_figures/NICD_tss_dist.png"
png(filename, width = 5, height = 5, units = "in", res = 300)
plotDistToTSS(NICD_peak_anno,
              title="Distribution of accessibility\nrelative to TSS")

dev.off()

print("DESeq2 Analysis and Peak Annotation done. Results are saved in the 'results/' directory.")



