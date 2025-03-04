# SET UP

library(DiffBind)
library(clusterProfiler)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(dplyr)
library(readr)


# READ IN DIFFBIND RESULTS: ----------------------------------------------------
diffbind_res_df <- read.csv(
  file.path("results/tables/diffbind", 
            paste0(project, "_diffbind_results.csv"))
  )

# Create dataframes for CTRL (Fold < 0) keeping only significant peaks 
# (FDR < 0.05)
# control_group name based on config_go_enrichment.R
assign(paste0(control_group, "_enrich"), diffbind_res_df %>% 
         dplyr::filter(FDR < 0.05 & Fold < 0) 
)

# Create dataframes for NICD keeping only significant peaks (FDR < 0.05)
# mutant_group name based on config_go_enrichment.R
assign(paste0(mutant_group, "_enrich"), diffbind_res_df %>% 
         dplyr::filter(FDR < 0.05 & Fold > 0)
)

# FUNCTIONAL ENRICHMENT ANALYSIS: ----------------------------------------------
# With clusterProfiler

# Make a list of genes for CTRL and NICD
#Extract the gene SYMBOL into a list

assign(paste0(control_group, "_genes"), 
       get(paste0(control_group, "_enrich"))$SYMBOL # get dataframe's SYMBOL col
)

assign(paste0(mutant_group, "_genes"), 
       get(paste0(mutant_group, "_enrich"))$SYMBOL # get dataframe's SYMBOL col
)


# GO term enrichment analysis for control_group
assign(paste0(control_group, "_GO_results"),
       enrichGO(gene = get(paste0(control_group, "_genes")), # get dataframe
                OrgDb = org.Mm.eg.db, # annotation database
                keyType = "SYMBOL", # gene ID type
                ont = "BP") # ontology: BP (biological process), MF (molecular function), CC (cellular component)
       )

# GO term enrichment analysis for mutant_group
assign(paste0(mutant_group, "_GO_results"),
       enrichGO(gene = get(paste0(mutant_group, "_genes")), # get dataframe
                OrgDb = org.Mm.eg.db, # annotation database
                keyType = "SYMBOL", # gene ID type
                ont = "BP") # ontology: BP (biological process), MF (molecular function), CC (cellular component)
)

# Write to csv files
filename <- paste0(control_group, "_GO_results_symbol.csv")

write.csv(
  get(paste0(control_group, "_GO_results")),   # specify dataframe, not string
  file = file.path(output_dir_tables, filename)
)

filename <- paste0(mutant_group, "_GO_results_symbol.csv")

write.csv(
  get(paste0(mutant_group, "_GO_results")),   # specify dataframe, not string
  file = file.path(output_dir_tables, filename)
)


print(paste0("GO results saved in ", output_dir_tables))





