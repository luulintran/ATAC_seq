# SET UP
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = FALSE)
BiocManager::install(c("DiffBind", 
                       "profileplyr", 
                       "ChIPseeker", 
                       "org.Mm.eg.db", 
                       "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                       "pheatmap",
                       "clusterProfiler",
                       "ReactomePA"), ask = FALSE)

library(DiffBind)
library(profileplyr)
library(clusterProfiler)
library(ChIPseeker)
library(ReactomePA)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(dplyr)
library(readr)

# Set working directory to project
# setwd("..")

# READ IN DATA:------------------------------------------------------------------------------------------
CTRL_enriched_list <- read.csv("results/tables/diffbind_tables/CTRL_e13_enriched_list.csv")
NICD_enriched_list <- read.csv("results/tables/diffbind_tables/NICD_e13_enriched_list.csv")

# FUNCTIONAL ENRICHMENT ANALYSIS: -----------------------------------------------------------------------------------------
# With clusterProfiler

# Make a list of genes for CTRL and NICD
#Extract the ENTREZ ID's into a list
CTRL_genes <- CTRL_enriched_list$geneId
NICD_genes <- NICD_enriched_list$geneId

# GO term enrichment analysis
GO_results_CTRL <- enrichGO(gene = CTRL_genes,
                            OrgDb = "org.Mm.eg.db", #annotation database
                            keyType = "ENTREZID", #gene id type
                            ont = "BP") #ontology: BP (biological process), MP (molecular function), CC (cellular component)

GO_results_NICD <- enrichGO(gene = NICD_genes,
                            OrgDb = "org.Mm.eg.db", #annotation database
                            keyType = "ENTREZID", #gene id type
                            ont = "BP") #ontology: BP (biological process), MF (molecular function), CC (cellular component)

# Save results as CSV files
write.csv(GO_results_CTRL, "results/tables/functional_enrichment_tables/GO_results_CTRL.csv")
write.csv(GO_results_NICD, "results/tables/functional_enrichment_tables/GO_results_NICD.csv")

# Make a barplot **********************************************************************************
filename = "results/figures/functional_enrichment_figures/GO_barplot_CTR.png"
png(filename, width = 10, height = 10, units = "in", res = 300)
CTRL_go_barplot <- plot(barplot(GO_results_CTRL, showCategory = 20))

dev.off()

filename = "results/figures/functional_enrichment_figures/GO_barplot_NICD.png"
png(filename, width = 10, height = 10, units = "in", res = 300)
NICD_go_barplot <- plot(barplot(GO_results_NICD, showCategory = 20))

dev.off()


# REACTOME PATHWAY ANALYSIS: ---------------------------------------------------------------------------------------------------------------
# CTRL
filename ="results/figures/functional_enrichment_figures/PA_dotplot_CTRL.png"
png(filename, width = 10, height = 10, units = "in", res = 300)
pathway_CTRL <- enrichPathway(CTRL_genes,
                              organism = "mouse")
dotplot(pathway_CTRL)

dev.off()

filename ="results/figures/functional_enrichment_figures/PA_dotplot_NICD.png"
png(filename, width = 10, height = 10, units = "in", res = 300)
pathway_NICD <- enrichPathway(NICD_genes,
                              organism = "mouse")
dotplot(pathway_NICD)

dev.off()

## KEGG PATHWAY ENRICHMENT ANALYSIS: ---------------------------------------------------------------------------------------------------------
filename = "results/figures/functional_enrichment_figures/KEGG_dotplot_CTRL.png"
png(filename, width = 10, height = 10, units = "in", res = 300)
pathway_KEGG_CTRL <- enrichKEGG(CTRL_genes,
                                organism = "mouse")
dotplot(pathway_KEGG_CTRL)

dev.off()

filename = "results/figures/functional_enrichment_figures/KEGG_dotplot_NICD.png"
png(filename, width = 10, height = 10, units = "in", res = 300)
pathway_KEGG_NICD <- enrichKEGG(NICD_genes,
                                organism = "mouse")
dotplot(pathway_KEGG_NICD)
dev.off()

print("Functional enrichment analysis complete. Results are in results/ directory")

