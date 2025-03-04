# Save GO results as dataframes
# ontrol_group
assign(paste0(control_group, "_GO_results_df"),
       as.data.frame(get(paste0(control_group, "_GO_results")))
)

# mutant_group
assign(paste0(mutant_group, "_GO_results_df"),
       as.data.frame(get(paste0(mutant_group, "_GO_results")))
)