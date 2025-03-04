library(clusterProfiler)

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
# With KEGG Enrichment Analysis

# Make a list of genes for CTRL and NICD
#Extract the gene ENTREZ ID's into a list (KEGG expects entrez not symbol)

assign(paste0(control_group, "_genes"), 
       get(paste0(control_group, "_enrich"))$geneId # get dataframe's SYMBOL col
)

assign(paste0(mutant_group, "_genes"), 
       get(paste0(mutant_group, "_enrich"))$geneId # get dataframe's SYMBOL col
)

## KEGG PATHWAY ENRICHMENT ANALYSIS: -------------------------------------------

# run Reactome Pathway Analysis
# control_group
assign(
  paste0(control_group, "_kegg"),  # assign results
  enrichKEGG(
    gene = as.character(get(paste0(control_group, "_genes"))), # extract gene list
    organism = "mouse" # organism
  )
)

# mutant_group
assign(
  paste0(mutant_group, "_kegg"),  # assign results
  enrichKEGG(
    gene = as.character(get(paste0(mutant_group, "_genes"))), # extract gene list
    organism = "mouse" # organism
  )
)

# SAVE RESULTS: ----------------------------------------------------------------
# Write to csv files
filename <- paste0(control_group, "_reactome_results.csv")

write.csv(
  get(paste0(control_group, "_kegg")),   # specify dataframe, not string
  file = file.path(output_dir_tables, filename)
)

filename <- paste0(mutant_group, "_kegg_results.csv")

write.csv(
  get(paste0(mutant_group, "_kegg")),   # specify dataframe, not string
  file = file.path(output_dir_tables, filename)
)


print(paste0("KEGG enrichment results saved in ", output_dir_tables))

