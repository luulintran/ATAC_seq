library(dplyr)
library(ggplot2)

# Save GO results as dataframes
# ontrol_group
assign(paste0(control_group, "_GO_results_df"),
       as.data.frame(get(paste0(control_group, "_GO_results")))
)

# mutant_group
assign(paste0(mutant_group, "_GO_results_df"),
       as.data.frame(get(paste0(mutant_group, "_GO_results")))
)

# GO TERMS LISTS: --------------------------------------------------------------
# Define the GO terms of interest
GO_terms_list_up <- c('gliogenesis', 
                    'glial cell differentiation', 
                    'myelination', 
                    'glial cell development', 
                    'oligodendrocyte development',
                    'smoothened signaling pathway',
                    'cilium assembly',
                    'cilium organization')

GO_terms_list_down <- c('regulation of neurogenesis', 
                     'dendrite development', 
                     'positive regulation of neurogenesis', 
                     'regulation of synapse structure or activity', 
                     'regulation of neuron differentiation')

# FILTER DATAFRAMES FOR GO TERMS: ----------------------------------------------
# filter dataframe for GO terms from list

# mutant_group
assign(
  paste0(mutant_group, "_GO_terms"),  
  get(paste0(mutant_group, "_GO_results_df")) %>%  # Get the dataframe
    dplyr::filter(Description %in% GO_terms_list_up)  # Filter based on GO terms
)

# control_group
assign(
  paste0(control_group, "_GO_terms"),  
  get(paste0(control_group, "_GO_results_df")) %>%  # Get the dataframe
    dplyr::filter(Description %in% GO_terms_list_down)  # Filter based on GO terms
)

# CALCULATE -LOG10(PADJUST): ---------------------------------------------------

# mutant_group
assign(
  paste0(mutant_group, "_GO_terms"),
  get(paste0(mutant_group, "_GO_terms")) %>%
    mutate(log_p.adjust = -log10(p.adjust))
)


# control_group
assign(
  paste0(control_group, "_GO_terms"),
  get(paste0(control_group, "_GO_terms")) %>%
    mutate(log_p.adjust = -log10(p.adjust))
)